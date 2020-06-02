#!/usr/bin/env python

from Bio.SeqIO import parse
import pandas as pd
from glob import glob
import os
from os.path import join as opj, basename, splitext
import numpy as np


# bin info

def contig_map(sm):
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
    cumulative = 0
    size = sum(lengths)
    for l in sorted(lengths):
        cumulative += l
        if float(cumulative) / size >= 0.5:
            return l


def bin_stats(f):
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
    return {'bp': size, 'GC': gc_f, 'Mbp': round(size_mb, 2), 'mean_contig': mean_l, 'median_contig': median_l,
            'min_contig': min_l, 'max_contig': max_l, 'n50': n50_l, 'contigs': len(contig_lengths)}


def calculate_bin_stats(files):
    stats = {}
    for f in files:
        name = basename(f)
        name, _ = splitext(name)
        stats[name] = bin_stats(f)
    return stats


def binning_stats(sm):
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

# bin annotation

def count_rrna(sm):
    df = pd.read_csv(sm.input, sep="\t", usecols=[0, 2, 8], header=None,
                     names=["contig", "type", "fields"])
    types = [x.split(";")[0].split("=")[-1] for x in df.fields]
    bins = [x.split(";")[-1].split("=")[-1] for x in df.fields]
    _df = pd.DataFrame(data={'rRNA_type': types, 'Bin_Id': bins})
    dfc = _df.reset_index().groupby(["Bin_Id", "rRNA_type"]).count()
    table = dfc.pivot_table(columns="rRNA_type", index="Bin_Id", fill_value=0)[
        "index"]
    table.index.name = table.columns.name = ""

    missing = set(["16S_rRNA", "23S_rRNA", "5S_rRNA"]).difference(table.columns)
    if len(missing) > 0:
        table = pd.merge(table,
                         pd.DataFrame(columns=missing, index=table.index, data=0),
                         left_index=True, right_index=True)
    table = table.loc[:, ["5S_rRNA", "16S_rRNA", "23S_rRNA"]]
    table.index.name = "Bin_Id"
    table.to_csv(sm.output, sep="\t", index=True)

def count_trna(sm):
    df = pd.read_csv(sm.input, sep="\t")
    dfc = df.groupby(["tRNA_type","Bin_Id"]).count().reset_index().loc[:,["tRNA_type","tRNA#","Bin_Id"]]
    table = dfc.pivot_table(columns="tRNA_type", index="Bin_Id")["tRNA#"].fillna(0)
    table.index.name = table.columns.name = ""
    table.to_csv(sm.output[0], sep="\t", index=True)

    total = {}
    for m in table.index:
        c = len(table.loc[m, table.loc[m]>0])
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


def generate_bin_list(input, outdir):
    genomes = []
    for f in input:
        if os.path.getsize(f) == 0:
            continue
        items = f.split("/")
        l, group, binner = items[-2], items[-3], items[-4]
        bindir = os.path.dirname(f)
        if binner == "concoct":
            bindir = opj(bindir, "fasta")
        # get absolute path for indir
        abs_in = os.path.abspath(bindir)
        df = pd.read_csv(f, sep="\t", index_col=0)
        uniq_suffix = "{group}.{l}".format(group=group, l=l)
        if not df.index.str.contains("metabat")[0]:
            df.rename(index = lambda x: "{binner}.{x}".format(binner=binner,
                                                              x=x), inplace=True)
        idmap = dict(zip(df.index, ["{x}.{s}".format(x=x, s=uniq_suffix) for x in df.index]))
        for bin_id, uniq_id in idmap.items():
            src = opj(abs_in, "{}.fa".format(bin_id))
            dst = opj(outdir, "{}.fa".format(uniq_id))
            if os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)
            genomes.append(dst)
    return genomes


def generate_ref_list(input, outdir):
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
    with open(output, 'w') as fh:
        for g in genomes:
            fh.write("{}\n".format(g))
    return genomes


def generate_fastANI_lists(sm):
    bins = generate_bin_list(sm.input.bins, sm.params.outdir)
    refs = generate_ref_list(sm.input.refs, sm.params.outdir)
    genomes = bins + refs
    write_list(genomes, sm.output[0])
    genomes.reverse()
    write_list(genomes, sm.output[1])


def fastani2dist(f, frac=0.5):
    """
    Reads the output.txt file from fastANI and generates a distance matrix

    :param f: Input file (out.txt typically)
    :param frac: Fraction of genome overlap to evaluate a pair of genomes
    :return: pandas DataFrame with distances
    """
    df = pd.read_csv(f, sep="\t", header=None,
                     names=["query", "ref", "ANI", "mapped_fragments",
                            "total_fragments"])
    # Add column with aligned fraction
    df = df.assign(frac=pd.Series(df.mapped_fragments / df.total_fragments,
                                  index=df.index))
    # Filter dataframe by aligned fraction
    df = df.loc[df.frac >= frac]
    # Rename genomes
    for key in ["query", "ref"]:
        df[key] = [x.split("/")[-1].replace(".fna", "").replace(".fa", "") for x
                   in df[key]]
    # Pivot
    table = df.pivot_table(index="query", columns="ref", values="ANI")
    # Calculate distance as 1-ANI/100
    dist = 1 - table.fillna(0).div(100)
    return dist


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
        edges = list(nx.dfs_edges(g,n))
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
        for j in range(i+1, len(dist_mat.columns)):
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
    dist = fastani2dist(sm.input[0], frac=sm.params.frac)
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
    toolbox = {"contig_map": contig_map,
               "count_tRNA": count_trna,
               "count_rRNA": count_rrna,
               "binning_stats": binning_stats,
               "download_ref_genome": download_ref_genome,
               "generate_fastANI_lists": generate_fastANI_lists,
               "cluster_genomes": cluster_genomes}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
