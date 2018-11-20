import pandas as pd, re, sys
from argparse import ArgumentParser

ec_re = re.compile("EC:[1-6\-].[0-9\-]+.[0-9\-]+.[0-9\-]+")


def read_ecmap(fh):
    enzymes = []
    proteins = []
    for line in fh:
        items = line.split("\t")
        m = ec_re.search(items[2])
        try:
            ec = m.group().split(":")[1]
        except AttributeError:
            continue
        member = "{}.{}".format(items[0], items[1])
        proteins.append(member)
        enzymes.append(ec)
    return enzymes, proteins


def ecmap(f):
    with open(f) as fh:
        enzymes, proteins = read_ecmap(fh)
    return enzymes, proteins


def read_cogmap(fh):
    cogs = []
    proteins = []
    for line in fh:
        items = line.split("\t")
        prots = items[-1].split(",")
        cog = [items[1]] * len(prots)
        cogs += cog
        proteins += prots
    return cogs, proteins


def cogmap(f):
    with open(f) as fh:
        cogs, proteins = read_cogmap(fh)
    return cogs, proteins


def cog2ec(map_df, frac=0.5):
    # Group by cog and enzyme to get number of each EC assignment per cog
    map_df_counts = map_df.groupby(["enzyme", "cog"]).count().reset_index()
    map_df_counts.index = map_df_counts.cog
    map_df_counts.drop("cog", axis=1, inplace=True)
    map_df_counts.sort_index(inplace=True)

    # Count total number of proteins per cog
    cog_counts = map_df_counts.groupby(level=0).sum()

    # Divide enzyme assignment number by total protein number to get fraction of each assignment
    ecfrac = map_df_counts.protein.div(cog_counts.protein).reset_index()
    # Get index of where fraction is above threshold
    index = ecfrac.loc[ecfrac.protein >= frac].index

    # Return mappings where fraction is above threshold
    return map_df_counts.iloc[index]


def main():
    parser = ArgumentParser()
    parser.add_argument("-c", "--conversion", type=str, required=True,
                        help="EGGNOG protein members to enzyme map file ('eggnog4.protein_id_conversion.tsv')")
    parser.add_argument("-m", "--members", type=str, required=True,
                        help="EGGNOG protein members file ('NOG.members.tsv')")

    args = parser.parse_args()

    enzymes, proteins = ecmap(args.conversion)
    ecmap_df = pd.DataFrame(data={"enzyme": enzymes, "protein": proteins})

    cogs, proteins = cogmap(args.members)
    cogmap_df = pd.DataFrame(data={"cog": cogs, "protein": proteins})

    map_df = pd.merge(ecmap_df, cogmap_df, left_on="protein", right_on="protein")

    cog2ec_df = cog2ec(map_df)

    cog2ec_df.loc[:, "enzyme"].to_csv(sys.stdout, sep="\t")


if __name__ == '__main__':
    main()