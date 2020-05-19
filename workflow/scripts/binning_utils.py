#!/usr/bin/env python

from Bio.SeqIO import parse
import pandas as pd
from glob import glob
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


def main(sm):
    toolbox = {"contig_map": contig_map,
               "count_tRNA": count_trna,
               "count_rRNA": count_rrna,
               "binning_stats": binning_stats}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
