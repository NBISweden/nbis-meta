#!/usr/bin/env python

from ete3 import NCBITaxa
import pandas as pd
from argparse import ArgumentParser


def get_seqids(df, seqid2taxid, ncbitaxa):
    filtered = pd.DataFrame()
    for i in df.index:
        taxid = df.loc[i,"taxID"]
        # Intersect descendants with those we have genomes for
        descendants = seqid2taxid.loc[seqid2taxid["taxID"].isin(ncbitaxa.get_descendant_taxa(taxid)+[taxid])]
        if len(descendants) == 0:
            continue
        filtered = pd.concat([filtered,descendants])
    return filtered


def filter_genomes(input,db,mapfile,min_read_count,output):
    df = pd.read_csv(input, sep="\t")
    # Filter to number of reads
    df = df.loc[df.numReads > min_read_count]
    # Filter to ranks genus, species and leaf
    df = df.loc[df.taxRank.isin(["leaf", "species", "genus"])]
    # Read the centrifuge sqlite database
    ncbitaxa = NCBITaxa(db)
    # Read the seqid -> taxid map file
    seqid2taxid = pd.read_csv(mapfile, header=None, names=["seq", "taxID"], index_col=0, sep="\t")
    filtered = get_seqids(df, seqid2taxid, ncbitaxa)
    if len(filtered) == 0:
        filtered = pd.DataFrame(columns=["taxID"], index=["seq"])
    filtered.to_csv(output, sep="\t", index=True)


def main(args):
    filter_genomes(args.input, args.db, args.mapfile, args.min_read_count, args.output)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-d", "--db", type=str, required=True)
    parser.add_argument("-m", "--min_read_count", type=int, required=True)
    parser.add_argument("-M", "--mapfile", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    args = parser.parse_args()
    main(args)
