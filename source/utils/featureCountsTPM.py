#!/usr/bin/env python

import pandas as pd
import sys
from argparse import ArgumentParser
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def calculate_tpm(df, readlength):
    sampleName = df.columns[-1]
    # 1. Calculate t for sample
    # t = (reads_mapped_to_gene * read_length) / length_of_gene
    # Multiply gene counts with read length, then divide by the Length column (in kb)
    logging.info("Normalizing by read length and gene length")
    t = df[sampleName].multiply(readlength).div(df["Length"].div(1000))
    df = df.assign(t=pd.Series(t, index=df.index))
    # 2. Calculate T
    # T = sum(t)
    logging.info("Calculating sum of normalized values")
    T = df["t"].sum()
    # 3. Calculate TPM
    # TPM = t*10^6 / T
    logging.info("Normalization factor T is {}".format(T))
    logging.info("Calculating TPM")
    TPM = (df["t"].multiply(1000000)).div(T)
    df = df.assign(TPM=pd.Series(TPM, index=df.index))
    return df


def main():
    p = ArgumentParser()
    p.add_argument("-i", "--infile", required=True,
                   help="Infile with counts from featureCounts")
    p.add_argument("--rl", dest="readlength", required=True, type=float,
                   help="Specify read length for sample being normalized")
    p.add_argument("--sampleName",
                   help="Specify sample name. If not specified, the basename will be extracted from the last column.")
    p.add_argument("-o", "--outfile",
                   help="Write TPM normalized results to outfile")
    p.add_argument("--rc", dest="rawcounts",
                   help="Also write the raw counts to this outfile")
    args = p.parse_args()
    logging.info("Reading counts from {}".format(args.infile))
    df = pd.read_table(args.infile, skiprows=1)
    if args.sampleName:
        sampleName = args.sampleName

    else:
        sampleName = df.columns[-1]
        sampleName = sampleName.split("/")[-1]

    if args.outfile:
        tpm_out = args.outfile
    else:
        tpm_out = sys.stdout

    df.columns = list(df.columns)[0:-1] + [sampleName]
    df = calculate_tpm(df, args.readlength)

    df_tpm = df.iloc[:, [0, -1]]
    df_tpm.columns = ["gene_id", sampleName]
    df_tpm.to_csv(tpm_out, sep="\t", index=False)

    if args.rawcounts:
        df_raw = df.iloc[:, [0, 6]]
        df_raw.columns = ["gene_id", sampleName]
        df_raw.to_csv(args.rawcounts, sep="\t", index=False)


if __name__ == '__main__':
    main()
