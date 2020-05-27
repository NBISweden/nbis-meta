#!/usr/bin/env python

import pandas as pd
from Bio.SeqIO import parse
from pathlib import Path


# file io

def fasta2bed(sm):
    with open(sm.input[0], 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for record in parse(fhin, "fasta"):
            fhout.write("{}\t{}\t{}\n".format(record.id, 0, len(record)))


# statistics

def store_lengths(f, minlen=False):
    """
    Reads lengths of contigs from fasta
    :param f: fasta file
    :param minlen: minimum length to store
    :return: pandas DataFrame of lengths
    """
    r = {}
    for record in parse(f, "fasta"):
        if minlen:
            if len(record.seq) < minlen:
                continue
        r[record.id] = len(record.seq)
    df = pd.DataFrame(r, index=["length"]).T
    return df


def size_distribute(df, lengths=None):
    """
    Calculates the distribution of an assembly in length bins

    For each <l> in <lengths> calculate for contigs >= <l>:
     n = the number of contigs
     s = the total length in bp
     p = the fraction of lengths / total assembly size
    :param df: pandas DataFrame of lengths
    :param lengths: intervals at which to calculate stats
    :return: pandas DataFrame
    """
    if lengths is None:
        lengths = [0, 100, 250, 500, 1000, 2500, 5000, 10000, 15000, 20000,
                   25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000,
                   125000, 150000, 200000, 250000, 500000]
    size_dist = {}
    for i, l in enumerate(lengths):
        if len(df.loc[df.length >= l]) == 0:
            break
        n = len(df.loc[df.length >= l])
        s = int(df.loc[df.length >= l].sum())
        p = int(df.loc[df.length >= l].sum()) / float(df.sum()) * 100
        size_dist[i] = {"min_length": l, "num_contigs": n, "total_length": s,
                        "%": p}
    size_dist_df = pd.DataFrame(size_dist).T
    size_dist_df = size_dist_df[
        ["min_length", "num_contigs", "total_length", "%"]]
    return size_dist_df


def calculate_n_stats(df):
    """
    Calculates n50 and n90 statistics from a list of lengths
    :param df: pandas DataFrame of contig lengths
    :return:
    """
    df.sort_values("length", inplace=True, ascending=True)
    size = int(df.sum())
    N50_length = N90_length = 0
    cumulative = 0
    for contig in df.index:
        l = df.loc[contig, "length"]
        cumulative += l
        if float(cumulative) >= 0.5 * size and not N50_length:
            N50_length = l
        elif float(cumulative) >= 0.1 * size and not N90_length:
            N90_length = l
    return N50_length, N90_length


def calculate_length_stats(df):
    """
    Calculates length statistics from a dataframe
    :param df: pandas DataFrame with contig lengths
    :return:
    """
    contigs = len(df)
    total_size = int(df.sum())
    min_length = int(df["length"].min())
    max_length = int(df["length"].max())
    avg_length = float(df["length"].mean())
    median_length = float(df["length"].median())
    return contigs, total_size, min_length, max_length, avg_length, median_length


def generate_stat_df(contig_lengths):
    """
    Generates statistics from a dataframe of contig lengths
    :param contig_lengths: pandas DataFrame
    :return:
    """
    index = ["contigs", "total_size_bp", "min_length", "max_length",
             "avg_length", "median_length", "N50_length", "N90_length"]
    stat_items = calculate_length_stats(contig_lengths)
    n50_length, n90_length = calculate_n_stats(contig_lengths)
    stat_df = pd.DataFrame([stat_items[0], stat_items[1], stat_items[2],
                            stat_items[3], stat_items[4], stat_items[5],
                            n50_length, n90_length], index=index).T
    return stat_df


def stats(sm):
    """
    Reads a list of assembly fasta files and generates statistics

    :param sm: snakemake object
    :return:
    """
    stat_result = pd.DataFrame()
    sizedist_result = pd.DataFrame()
    for f in sm.input.fa:
        p = Path(f)
        name = p.parent.name
        contig_lengths = store_lengths(f)
        stat_df = generate_stat_df(contig_lengths)
        size_dist = size_distribute(contig_lengths)
        stat_df["assembly"] = [name]*len(stat_df)
        size_dist["assembly"] = [name]*len(size_dist)
        stat_result = pd.concat([stat_result, stat_df])
        sizedist_result = pd.concat([sizedist_result,size_dist])
    stat_result = stat_result[["assembly", "contigs", "total_size_bp",
                               "min_length", "max_length", "avg_length",
                               "median_length", "N50_length", "N90_length"]]
    stat_result.to_csv(sm.output[0], sep="\t", index=False)
    sizedist_result = sizedist_result[["assembly", "min_length",
                                       "num_contigs", "total_length", "%"]]
    sizedist_result.to_csv(sm.output[1], sep="\t", index=False)

# assembly input

def metaspades_input(sm):
    """
    Generates fastq files to use as input for metaspades assembler

    :param sm: snakemake object
    :return:
    """
    from common import rename_records
    files = {"R1": [], "R2": [], "se": []}
    assembly_dict = sm.params.assembly
    # Collect all files belonging to the assembly group
    for sample in assembly_dict.keys():
        for unit in assembly_dict[sample]:
            for pair in assembly_dict[sample][unit].keys():
                files[pair].append(
                    assembly_dict[sample][unit][pair][0])
    # Rename and concatenate reads (required for Metaspades)
    with open(sm.output.R1, 'w') as fh1, open(sm.output.R2, 'w') as fh2, open(
        sm.output.se, 'w') as fhse:
        i = 0
        for f in files["R1"]:
            f2 = files["R2"][i]
            fh1 = rename_records(f, fh1, i)
            fh2 = rename_records(f2, fh2, i)
            i += 1
        for i, f in enumerate(files["se"], start=i):
            fhse = rename_records(f, fhse, i)


def megahit_input(sm):
    """
    Genereate input lists for megahit assembler

    :param sm: snakemake object
    :return:
    """
    files = {"R1": [], "R2": [], "se": []}
    assembly_dict = sm.params.assembly
    for sample in assembly_dict.keys():
        for unit in assembly_dict[sample]:
            for pair in assembly_dict[sample][unit].keys():
                files[pair].append(assembly_dict[sample][unit][pair][0])
    with open(sm.output.R1, 'w') as fh1, \
        open(sm.output.R2, 'w') as fh2, \
        open(sm.output.se, 'w') as fhse:
        fh1.write(",".join(files["R1"]))
        fh2.write(",".join(files["R2"]))
        fhse.write(",".join(files["se"]))


def main(sm):
    toolbox = {"assembly_stats": stats,
               "generate_megahit_input": megahit_input,
               "generate_metaspades_input": metaspades_input,
               "fasta2bed": fasta2bed}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
