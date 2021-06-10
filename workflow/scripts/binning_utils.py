#!/usr/bin/env python

import pandas as pd
from glob import glob
import os
from os.path import join as opj, basename, splitext
import numpy as np


# bin info

def contig_map(sm):
    """
    Generates a map of bin->contig id
    :param sm: snakemake object
    :return:
    """
    from Bio.SeqIO import parse
    files = glob(opj(sm.params.dir, "*.fa"))
    if len(files) == 0:
        with open(sm.output[0], 'w') as fhout:
            pass
        return
    for f in files:
        with open(f, 'r') as fhin, open(sm.output[0], 'w') as fhout:
            genome, _ = splitext(basename(f))
            for record in parse(fhin, "fasta"):
                fhout.write("{}\t{}\n".format(genome, record.id))


# bin stats

def n50(lengths):
    """
    Calculate N50 stats
    :param lengths:
    :return:
    """
    cumulative = 0
    size = sum(lengths)
    for l in sorted(lengths):
        cumulative += l
        if float(cumulative) / size >= 0.5:
            return l


def bin_stats(f):
    """
    Generate bin statistics
    :param f:
    :return:
    """
    from Bio.SeqIO import parse
    size = 0
    gc = 0
    contig_lengths = []
    for record in parse(f, "fasta"):
        l = len(record)
        seq = record.seq.upper()
        g = seq.count("G")
        c = seq.count("C")
        gc += g + c
        size += l
        contig_lengths.append(l)
    gc_f = round(float(gc) / size * 100, 2)
    size_mb = size / 1000000
    mean_l = round(np.mean(contig_lengths), 2)
    median_l = round(np.median(contig_lengths), 2)
    min_l = np.min(contig_lengths)
    max_l = np.max(contig_lengths)
    n50_l = n50(contig_lengths)
    return {'bp': size, 'GC': gc_f, 'Mbp': round(size_mb, 2),
            'mean_contig': mean_l, 'median_contig': median_l,
            'min_contig': min_l, 'max_contig': max_l, 'n50': n50_l,
            'contigs': len(contig_lengths)}


def calculate_bin_stats(files):
    """
    Calls bin statistics for each file
    :param files:
    :return:
    """
    stats = {}
    for f in files:
        name = basename(f)
        name, _ = splitext(name)
        stats[name] = bin_stats(f)
    return stats


def binning_stats(sm):
    """
    Main function for calculating bin statistics
    :param sm:
    :return:
    """
    files = glob(opj(sm.params.dir, "*.fa"))
    if len(files) == 0:
        with open(sm.output[0], 'w') as fh:
            fh.write("No bins found\n")
        return
    else:
        stats = calculate_bin_stats(files)
        cols = ["bp", "Mbp", "GC", "contigs", "n50", "mean_contig",
                "median_contig", "min_contig", "max_contig"]
        df = pd.DataFrame(stats).T[cols]
        df.index.name = "bin"
        df.sort_values("bp", ascending=False, inplace=True)
        df.to_csv(sm.output[0], sep="\t", index=True, header=True)


# checkm utils

def remove_checkm_zerocols(sm):
    """
    Reads checkm coverage and removes samples with no reads mapped to contigs

    :param sm:
    :return:
    """
    df = pd.read_csv(sm.input[0], header=0, sep="\t")
    if df.shape[0] == 0:
        with open(sm.output[0], 'w') as fh:
            fh.write("NO BINS FOUND\n")
        return
    # base columns are independent of samples
    base_columns = ["Sequence Id", "Bin Id", "Sequence length (bp)"]
    # get all columns with mapped read counts
    cols = [x for x in df.columns if "Mapped reads" in x]
    # sum counts
    df_sum = df.loc[:, cols].sum()
    # get columns with zero reads mapped
    zero_cols = df_sum.loc[df_sum == 0].index
    cols_to_drop = []
    # find suffix of zero cols
    for c in list(zero_cols):
        suffix = ".{}".format(c.split(".")[-1])
        if suffix == ".Mapped reads":
            suffix = ""
        cols_to_drop += (
            "Mapped reads{s},Bam Id{s},Coverage{s}".format(s=suffix)).split(",")
    df_checked = df.drop(cols_to_drop, axis=1)
    # check that there are remaining columns
    diff_cols = set(df_checked.columns).difference(base_columns)
    if len(diff_cols) > 0:
        df_checked.to_csv(sm.output[0], sep="\t", index=False, header=True)
    else:
        with open(sm.output[0], 'w') as fhout:
            pass

# bin annotation

def count_rrna(sm):
    """
    Counts rRNA genes in bins
    :param sm:
    :return:
    """
    df = pd.read_csv(sm.input[0], sep="\t", usecols=[0, 2, 8], header=None,
                     names=["contig", "type", "fields"])
    # If empty dataframe, just write an empty file
    if df.shape[0] == 0:
        table = pd.DataFrame(columns=["5S_rRNA", "16S_rRNA", "23S_rRNA"])
        table.index.name = "Bin Id"
        table.to_csv(sm.output[0], sep="\t", header=True)
        return
    types = [x.split(";")[0].split("=")[-1] for x in df.fields]
    bins = [x.split(";")[-1].split("=")[-1] for x in df.fields]
    _df = pd.DataFrame(data={'rRNA_type': types, 'Bin_Id': bins})
    dfc = _df.reset_index().groupby(["Bin_Id", "rRNA_type"]).count()
    table = dfc.pivot_table(columns="rRNA_type", index="Bin_Id", fill_value=0)[
        "index"]
    table.index.name = table.columns.name = ""

    missing = set(["16S_rRNA", "23S_rRNA", "5S_rRNA"]).difference(table.columns)
    if len(missing) > 0:
        table = pd.merge(table, pd.DataFrame(columns=missing, index=table.index,
                                             data=0), left_index=True,
                         right_index=True)
    table = table.loc[:, ["5S_rRNA", "16S_rRNA", "23S_rRNA"]]
    table.index.name = "Bin_Id"
    table.to_csv(sm.output[0], sep="\t", index=True)


def count_trna(sm):
    """
    Counts tRNA genes in bins
    :param sm:
    :return:
    """
    df = pd.read_csv(sm.input[0], sep="\t")
    dfc = df.groupby(["tRNA_type", "Bin_Id"]).count().reset_index().loc[:,
          ["tRNA_type", "tRNA#", "Bin_Id"]]
    table = dfc.pivot_table(columns="tRNA_type", index="Bin_Id")[
        "tRNA#"].fillna(0)
    table.index.name = table.columns.name = ""
    table.to_csv(sm.output[0], sep="\t", index=True)

    total = {}
    for m in table.index:
        c = len(table.loc[m, table.loc[m] > 0])
        total[m] = c
    table = pd.DataFrame(total, index=["tRNAs"]).T
    table.index.name = "Bin_Id"
    table.to_csv(sm.output[1], sep="\t", index=True)


# genome clustering

def fetch_genome(ftp_base, outfile):
    """
    Performs the ftp fetching
    :param ftp_base: Base url to RefSeq/GenBank ftp for the genome
    :param outfile: Output filename
    :return:
    """
    from urllib import request
    ftp_base = ftp_base.rstrip("/")
    n = os.path.basename(ftp_base)
    fna = opj(ftp_base, "{}_genomic.fna.gz".format(n))
    r = request.urlretrieve(fna, outfile)
    return r


def unzip(z, o):
    """
    Unzips a zipped fasta file and removes the zipped file
    :param z: input gzip file
    :param o: output plain text file
    :return:
    """
    import gzip as gz
    with gz.open(z, 'rt') as fhin, open(o, 'w') as fhout:
        fhout.write(fhin.read())
    os.remove(z)


def download_ref_genome(sm):
    """
    Fetches reference genomes for clustering with fastANI
    :param sm: snakemake object
    :return:
    """
    fetch_genome(sm.params.ftp_base, "{}.gz".format(sm.output[0]))
    unzip("{}.gz".format(sm.output[0]), sm.output[0])


def generate_bin_list(input, outdir, min_completeness, max_contamination):
    """
    Generates a list of bins to use for fastANI
    Also symlinks each bin file into the fastANI folder
    :param input: List of checkm genome statistics files
    :param outdir: Output directory to store symlinks
    :param min_completeness: Minimum completeness level (%)
    :param max_contamination: Maximum contamination level (%)
    :return:
    """
    genomes = []
    for f in input:
        with open(f, 'r') as fh:
            if fh.readline().rstrip() == "NO BINS FOUND":
                continue
        items = f.split("/")
        # extract wildcards from file path
        binner, assembly, l = items[-5], items[-4], items[-3]
        bindir = os.path.dirname(os.path.dirname(f))
        # get absolute path for bin directory
        abs_in = os.path.abspath(bindir)
        # read the checkm summary file
        df = pd.read_csv(f, sep="\t", index_col=0)
        # filter to at least <min_completeness> and at most <max_contamination>
        df = df.loc[(df.Completeness >= min_completeness) & (df.Contamination <= max_contamination)]
        if df.shape[0] == 0:
            continue
        # generate a unique suffix for each bin
        uniq_suffix = "{assembly}.{l}".format(assembly=assembly, l=l)
        # make a map of the bin id and the unique suffix
        idmap = dict(zip(df.index,
                         ["{x}.{s}".format(x=x, s=uniq_suffix) for x in
                          df.index]))
        # create symlink in the output path for each bin id that points
        # to the original fasta file
        for bin_id, uniq_id in idmap.items():
            src = opj(abs_in, "{}.fa".format(bin_id))
            dst = opj(outdir, "{}.fa".format(uniq_id))
            if os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)
            genomes.append(dst)
    return genomes


def generate_ref_list(input, outdir):
    """
    Generates a list of reference genomes
    Also symlinks each reference file into the fastANI folder

    :param input:
    :param outdir:
    :return:
    """
    genomes = []
    for f in input:
        basename = os.path.basename(f)
        src = os.path.abspath(f)
        dst = opj(outdir, basename)
        if os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)
        genomes.append(dst)
    return genomes


def write_list(genomes, output):
    """
    Writes the fastANI lists to file
    :param genomes: list of genomes to use for fastANI
    :param output:
    :return:
    """
    with open(output, 'w') as fh:
        for g in genomes:
            fh.write("{}\n".format(g))
    return genomes


def generate_fastANI_lists(sm):
    """
    Main function for generating fastANI lists
    :param sm:
    :return:
    """
    bins = generate_bin_list(sm.input.bins, sm.params.outdir,
                             sm.params.completeness,
                             sm.params.contamination)
    refs = generate_ref_list(sm.input.refs, sm.params.outdir)
    genomes = bins + refs
    write_list(genomes, sm.output[0])
    genomes.reverse()
    write_list(genomes, sm.output[1])


def check_pairs(pairs, min_frags):
    allowed_pairs = {}
    for i in pairs.index:
        q = pairs.loc[i, "query"]
        r = pairs.loc[i, "ref"]
        for key in [q, r]:
            if not key in allowed_pairs.keys():
                allowed_pairs[key] = []
        if pairs.loc[i, "aligned"] >= min_frags:
            allowed_pairs[q].append(r)
            allowed_pairs[r].append(q)
    return allowed_pairs


def fastani2dist(mat, txt, min_frags):
    """
    Converts the fastANI out.txt.matrix file to a pandas DataFrame
    :param mat: Distance matrix file
    :param txt: Pairwise output table with ANI and aligned + total fragments
    :param min_frags: Minimum aligned fragments to compare two genomes
    :return:
    """
    # read the pairwise table
    pairs = pd.read_csv(txt, header=None,
                        sep="\t",
                        names=["query", "ref", "ANI", "aligned", "total"])
    for key in ["query", "ref"]:
        pairs[key] = [x.split("/")[-1].replace(".fna", "").replace(".fa", "")
                      for x in pairs[key]]
    allowed_pairs = check_pairs(pairs, min_frags)
    genomes = list(pd.read_table(mat, index_col=0, skiprows=1, sep="\t",
                                 header=None, usecols=[0]).index)
    genomes = [os.path.splitext(os.path.basename(g))[0] for g in genomes]
    r = {}
    with open(mat, 'r') as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            line = line.rstrip()
            items = line.rsplit("\t")
            genome = os.path.splitext(os.path.basename(items[0]))[0]
            r[genome] = {}
            for j, item in enumerate(items[1:]):
                genome2 = genomes[j]
                # check that the pairing is allowed
                if item == "NA" or genome not in allowed_pairs[genome2] or genome2 not in allowed_pairs[genome]:
                    item = np.nan
                else:
                    item = float(item)
                r[genome][genome2] = item
    df = pd.DataFrame(r)
    df.fillna(0, inplace=True)
    return 1-df.div(100)


def cluster(linkage):
    """
    Cluster all genomes based on established linkages using networkx

    :param linkage: Dictionary of dictionaries
    :return: A dictionary with cluster index and list of genomes
    """
    import networkx as nx
    g = nx.from_dict_of_dicts(linkage)
    clustered = []
    clusters = {}
    clust_num = 1
    for n in g.nodes():
        c = [n]
        if n in clustered: continue
        edges = list(nx.dfs_edges(g, n))
        for e in edges:
            n1, n2 = e
            clustered += [n1, n2]
            c += [n1, n2]
        c = list(set(c))
        clusters[clust_num] = c[:]
        clust_num += 1
    return clusters


def generate_linkage(dist_mat, max_dist):
    """
    Create a nested dictionary linking genomes if their distance is within
    a certain threshold.

    :param dist_mat: pandas DataFrame with distances
    :param max_dist: maximum allowed distance to link genomes
    :return: a nested dictionary
    """
    linkage = {}
    for i in range(len(dist_mat.index)):
        g1 = dist_mat.index[i]
        if not g1 in linkage.keys():
            linkage[g1] = {}
        for j in range(i + 1, len(dist_mat.columns)):
            g2 = dist_mat.columns[j]
            if not g2 in linkage.keys():
                linkage[g2] = {}
            distance = dist_mat.iloc[i, j]
            if distance <= max_dist:
                linkage[g1][g2] = ""
                linkage[g2][g1] = ""
    return linkage


def cluster_genomes(sm):
    """
    Main function to run clustering of genomes

    :param sm: snakemake object
    :return:
    """
    dist = fastani2dist(sm.input.mat, sm.input.txt, sm.params.minfrags)
    linkage = generate_linkage(dist, sm.params.thresh)
    clusters = cluster(linkage)
    write_clusters(clusters, sm.output[0])


def write_clusters(clusters, outfile):
    """
    Sorts clusters by size and writes to file
    :param clusters: Dictionary of clusters with genomes as a list
    :param outfile: Output file path
    :return:
    """
    import operator
    # Calculate cluster sizes
    cluster_sizes = {}
    for clust_num, l in clusters.items():
        cluster_sizes[clust_num] = len(l)
    # Sort clusters by sizes
    sorted_clusters = sorted(cluster_sizes.items(), key=operator.itemgetter(1),
                             reverse=True)
    # Write table
    with open(outfile, 'w') as fh:
        for i, item in enumerate(sorted_clusters, start=1):
            old_num = item[0]
            for g in clusters[old_num]:
                fh.write("Cluster{}\t{}\n".format(i, g))


def main(sm):
    toolbox = {"contig_map": contig_map, "count_tRNA": count_trna,
               "remove_checkm_zerocols": remove_checkm_zerocols,
               "count_rRNA": count_rrna, "binning_stats": binning_stats,
               "download_ref_genome": download_ref_genome,
               "generate_fastANI_lists": generate_fastANI_lists,
               "cluster_genomes": cluster_genomes}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
