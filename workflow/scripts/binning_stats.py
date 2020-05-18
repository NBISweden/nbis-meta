#!/usr/bin/env python

from Bio.SeqIO import parse
import pandas as pd
from glob import glob
from os.path import join as opj
from os.path import basename
import numpy as np


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


def calculate_bin_stats(files, suffix):
    stats = {}
    for f in files:
        name = basename(f)
        name = name.rstrip(suffix)
        stats[name] = bin_stats(f)
    return stats


def main(sm):
    files = glob(opj(sm.params.dir, "*{}".format(sm.params.suffix)))
    if len(files) == 0:
        with open(sm.output[0], 'w') as fh:
            fh.write("No bins found\n")
        return
    else:
        stats = calculate_bin_stats(files, sm.params.suffix)
        cols = ["bp", "Mbp", "GC", "contigs", "n50", "mean_contig",
                "median_contig", "min_contig", "max_contig"]
        df = pd.DataFrame(stats).T[cols]
        df.index.name = "bin"
        df.sort_values("bp", ascending=False, inplace=True)
        df.to_csv(sm.output[0], sep="\t", index=True, header=True)


if __name__ == '__main__':
    main(snakemake)
