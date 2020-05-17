#!/usr/bin/env python

import pandas as pd
from Bio.SeqIO import parse
from pathlib import Path


def store_lengths(f, minlen=False):
    r = {}
    for record in parse(f, "fasta"):
        if minlen:
            if len(record.seq) < minlen:
                continue
        r[record.id] = len(record.seq)
    df = pd.DataFrame(r, index=["length"]).T
    return df


def size_distribute(df,
                    lengths=[0, 100, 250, 500, 1000, 2500, 5000, 10000, 15000,
                             20000, 25000, 30000, 35000, 40000, 45000, 50000,
                             75000, 100000, 125000, 150000, 200000, 250000,
                             500000]):
    size_dist = {}
    for i, l in enumerate(lengths):
        if len(df.loc[df.length >= l]) == 0:
            break
        n, s, p = len(df.loc[df.length >= l]), int(
            df.loc[df.length >= l].sum()), int(
            df.loc[df.length >= l].sum()) / float(df.sum()) * 100
        size_dist[i] = {"min_length": l, "num_contigs": n, "total_length": s,
                        "%": p}
    size_dist_df = pd.DataFrame(size_dist).T
    size_dist_df = size_dist_df[
        ["min_length", "num_contigs", "total_length", "%"]]
    return size_dist_df


def calculate_N_stats(df):
    df.sort_values("length", inplace=True, ascending=True)
    size = int(df.sum())
    N50 = N90 = False
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
    contigs = len(df)
    total_size = int(df.sum())
    min_length = int(df["length"].min())
    max_length = int(df["length"].max())
    avg_length = float(df["length"].mean())
    median_length = float(df["length"].median())
    return contigs, total_size, min_length, max_length, avg_length, median_length


def generate_stat_df(contig_lengths):
    index = ["contigs", "total_size_bp", "min_length", "max_length",
             "avg_length", "median_length", "N50_length", "N90_length"]
    contigs, total_size, min_length, max_length, avg_length, median_length = calculate_length_stats(
        contig_lengths)
    N50_length, N90_length = calculate_N_stats(contig_lengths)

    stat_df = pd.DataFrame([contigs, total_size, min_length, max_length,
                            avg_length, median_length, N50_length, N90_length],
                           index=index).T
    return stat_df

stat_result = pd.DataFrame()
sizedist_result = pd.DataFrame()
for f in snakemake.input:
    p = Path(f)
    name = p.parent.name
    contig_lengths = store_lengths(f)
    stat_df = generate_stat_df(contig_lengths)
    size_dist = size_distribute(contig_lengths)
    stat_df["assembly"] = [name]*len(stat_df)
    size_dist["assembly"] = [name]*len(size_dist)
    stat_result = pd.concat([stat_result,stat_df])
    sizedist_result = pd.concat([sizedist_result,size_dist])
stat_result = stat_result[["assembly","contigs", "total_size_bp",
                           "min_length", "max_length", "avg_length",
                           "median_length","N50_length", "N90_length"]]
stat_result.to_csv(snakemake.output[0], sep="\t", index=False)
sizedist_result = sizedist_result[["assembly", "min_length",
                                   "num_contigs", "total_length", "%"]]
sizedist_result.to_csv(snakemake.output[1], sep="\t", index=False)