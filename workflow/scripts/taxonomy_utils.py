#!/usr/bin/env python

import pandas as pd


def add_lower(d, ranks):
    """
    Propagates assignments from higher to lower taxonomic ranks,
    and adds a 'Unclassified.' prefix.
    :param d: dictionary of taxonomic assignments
    :param ranks: ranks for which to propagate
    :return:
    """
    last_known = d[ranks[0]]
    for rank in ranks[1:]:
        if d[rank] != "unknown":
            last_known = d[rank]
        else:
            if last_known == "Unclassified":
                d[rank] = last_known
            else:
                d[rank] = f"Unclassified.{last_known}"
    return d


def extract_lineage(d, ranks):
    taxonomy = {}
    for contig, lineage in d.items():
        taxonomy[contig] = {}
        for i, taxname in enumerate(lineage.split(";")):
            taxonomy[contig][ranks[i]] = taxname
        taxonomy[contig] = add_lower(taxonomy[contig], ranks)
    return taxonomy


def assign_orfs(sm):
    """
    Transfers taxonomic assignments from contigs down to ORFs called on contigs
    :param sm: snakemake object
    :return:
    """
    ranks = sm.params.ranks
    gff_df = pd.read_csv(
        sm.input.gff,
        header=None,
        sep="\t",
        comment="#",
        usecols=[0, 8],
        names=["contig", "id"],
    )
    # Extract ids
    ids = [
        "{}_{}".format(
            gff_df.loc[i, "contig"], gff_df.loc[i, "id"].split(";")[0].split("_")[-1]
        )
        for i in gff_df.index
    ]
    gff_df.loc[:, "id"] = ids
    # Read taxonomy for contigs
    tax_df = pd.read_csv(
        sm.input.tax,
        header=None,
        index_col=0,
        sep="\t",
        names=[
            "contig",
            "taxid",
            "rank",
            "assignment",
            "n_frags",
            "label_frags",
            "agree_frags",
            "support",
            "lineage",
            "short_rank_lineage",
        ],
    )
    # keep only contigs with assignments
    tax_df = tax_df.loc[tax_df.lineage == tax_df.lineage]
    taxdict = extract_lineage(tax_df.loc[:, "lineage"].to_dict(), ranks)
    taxdf = pd.DataFrame(taxdict).T
    # Merge dataframes
    orf_tax_df = pd.merge(
        gff_df, taxdf, left_on="contig", right_index=True, how="outer"
    )
    # When using 'outer' merging there may be contigs with no called ORF
    # but with a tax assignment. Drop these contigs.
    orf_tax_df = orf_tax_df.loc[orf_tax_df["id"] == orf_tax_df["id"]]
    # Set Unclassified for NA values
    orf_tax_df.fillna("Unclassified", inplace=True)
    # Set index to ORF ids
    orf_tax_df.set_index("id", inplace=True)
    orf_tax_df.drop("contig", axis=1, inplace=True)
    orf_tax_df.to_csv(sm.output.tax, sep="\t", index=True, header=True)


def main(sm):
    toolbox = {
        "assign_orfs": assign_orfs,
    }
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
