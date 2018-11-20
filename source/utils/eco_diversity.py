import skbio, os, numpy as np, pandas as pd, sys
from ete3 import NCBITaxa
from skbio.tree import TreeNode
from argparse import ArgumentParser

def read_taxonomy_counts(sample_path, rank):
    counts = {}
    with open(sample_path, 'r') as f:
        for l in f:
            l = l.strip('\n').split('\t')
            if l[3] == rank and not l[-1].lstrip() == 'unclassified':
                counts[l[-1].lstrip()] = int(l[1])
    return counts


def form_dataframe(sample_filenames, rank):
    """Form a pandas dataframe containing the counts.
    This can also be useful later for additional plots.
    Example:
                                     DA134_001  DA134_002  DA139_002
    Acanthamoeba castellanii           1.0        2.0        2.0
    Acaricomes phytoseiuli             8.0        6.0        NaN
    Acetobacter nitrogenifigens        4.0        NaN        NaN
    Acetobacter papayae                1.0        1.0        NaN
    Acetobacter pomorum                NaN        NaN        2.0
    """
    counts = {}
    for sample_path in sample_filenames:
        sample_name = os.path.basename(sample_path)
        counts[sample_name] = read_taxonomy_counts(sample_path, rank)
    df = pd.DataFrame.from_dict(counts, orient='columns', dtype=np.int64)
    df = df.fillna(value=np.int64(0))
    return df


def compute_alpha_diversities(df, ml):
    """
    Compute the alpha diversities.
    :param df: input count dataframe, as generater with form_dataframe()
    :param ml: measures list, should come from .rules
    :returns: a pandas dataframe with the different alpha diversities
    """
    # The numpy matrix oc counts, rows are sample counts
    mt = df.values.T.astype(np.dtype('int64'))
    # Sample list
    sl = list(df.columns.values)
    divs = {}
    for measure in ml:
        divs[measure] = skbio.diversity.alpha_diversity(measure, mt, ids=sl)
    df_adiv = pd.DataFrame.from_dict(divs, orient='columns')
    return df_adiv


def load_taxonomy_tree(otu_list):
    ncbi = NCBITaxa()
    # downloads the NCBI locally
    # TODO: it should maybe be done as a separate rule
    # ncbi.update_taxonomy_database()

    sp2taxid = ncbi.get_name_translator(otu_list)
    # The ETE docs scared me into checking if there are no or different matches
    for sp in sp2taxid:
        if len(sp2taxid[sp]) > 1:
            print("More than two NCBI taxonomy matches for:\n",
                  sp, sp2taxid[sp])
            sys.exit()
        if len(sp2taxid[sp]) < 1:
            print("No NCBI taxonomy match for:\n", sp, sp2taxid[sp])
            sys.exit()

    lineages = {}
    i = 1
    ncbi = NCBITaxa()
    for sp in sp2taxid:
        lineage = ncbi.get_lineage(sp2taxid[sp][0])
        names = ncbi.get_taxid_translator(lineage)
        lineages[sp] = [names[taxid] for taxid in lineage]
        i += 1
    tree = TreeNode.from_taxonomy(lineages.items())

    # Branches are required to have a length
    # TODO: Length should be parametrized
    for node in tree.postorder():
        node.length = 1

    # Lineages for species that are not found are silently not reported
    # They break the unifraq and have to be discarded from the dataframe
    # TODO: The notfound list should be reported!
    tips = set([tip.name for tip in tree.tips()])
    notfound = set([otu for otu in otu_list if otu not in tips])

    return tree, notfound


def compute_beta_unifraq(df, m):
    from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac

    # get the phylogenetic tree and drop OTUs that are not NCBI annotated
    otu_ids = list(df.index.values)
    tree, notfound = load_taxonomy_tree(otu_ids)
    df = df.drop(list(notfound))

    # The numpy matrix of counts, in which the rows are sample counts
    mt = df.values.T.astype(np.dtype('int64'))
    # Sample list
    sl = list(df.columns.values)
    # OTU list
    otu_ids = list(df.index.values)

    # The beta diversity matrix
    nsamples = len(sl)
    bm = np.zeros((nsamples, nsamples))
    # Compute the pairwise unifraq
    for i in range(nsamples):
        for j in range(i):
            u_counts = mt[i]
            v_counts = mt[j]
            if m == "unifraq":
                uu = unweighted_unifrac(u_counts, v_counts, otu_ids, tree)
            if m == "wunifraq":
                uu = weighted_unifrac(u_counts, v_counts, otu_ids, tree)
            bm[i, j] = uu
            bm[j, i] = uu
    return bm


def compute_beta_braycurtis(df):
    from skbio.diversity import beta_diversity
    l = []
    ids = []
    for i in range(0, len(df.columns)):
        l.append(df.iloc[:, i])
        ids.append(df.columns[i])
    array = np.array(l)
    counts = array.astype(int)
    r = beta_diversity(metric="braycurtis", counts = counts, ids = ids)
    bm = []
    r_df = r.to_data_frame()
    for i in range(0, len(r_df)):
        bm.append(list(r_df.iloc[i, :]))
    return bm

def compute_beta_diversities(df, ml):
    """
    Compute the beta diversities.
    :param df: counts dataframe, as generated with form_dataframe()
    :param ml: measures list, should come from .rules
    :returns: a pandas dataframe with the different beta diversities
    """
    divs = {}
    # Sample list
    sl = list(df.columns.values)
    # each divs[measure] is a numpy 2d array so I cannot make a dataframe
    # thus it is better to save it directly
    for measure in ml:
        if measure == "unifraq":
            divs[measure] = compute_beta_unifraq(df, measure)
        if measure == "wunifraq":
            divs[measure] = compute_beta_unifraq(df, measure)
        if measure == "braycurtis":
            divs[measure] = compute_beta_braycurtis(df)
    return divs, sl

def save_alpha_diversity(adiv, fh):
    adiv.to_csv(fh)
    return


def save_beta_diversity(bdiv, sl, fh):
    fh.write('Metric\t'+'\t'.join(sl)+'\n')
    for metric in bdiv:
        #fh.write('#' + metric + ':\n')
        for row in bdiv[metric]:
            fh.write(metric+"\t"+'\t'.join([str(i) for i in row])+'\n')
    return


def main():

    parser = ArgumentParser()
    parser.add_argument("-i", "--input", nargs = "+", required = True,
                        help="Input report files from classifier.")
    parser.add_argument("-r", "--rank", default="S",
                        help="Taxonomic rank at which to calculate diversity. Defaults to 'S' (species).")
    parser.add_argument("--div_type", type=str, required=True,
                        help="Type of diversity to calculate ('alpha' or 'beta')")
    parser.add_argument("-b", "--beta_measures", nargs="+", default=['braycurtis', 'wunifraq', 'unifraq'],
                        help="Type of diversity measures to calculate for beta-diversity")
    parser.add_argument("-a", "--alpha_measures", nargs="+", default=['shannon', 'simpson', 'observed_otus'],
                        help="Type of diversity measures to calculate for alpha-diversity")
    args = parser.parse_args()

    df = form_dataframe(args.input, args.rank)

    if args.div_type == "alpha":
        adiv = compute_alpha_diversities(df, args.alpha_measures)
        save_alpha_diversity(adiv, sys.stdout)

    elif args.div_type == "beta":
        bdiv, sl = compute_beta_diversities(df, args.beta_measures)
        save_beta_diversity(bdiv, sl, sys.stdout)

if __name__ == '__main__':
    main()
