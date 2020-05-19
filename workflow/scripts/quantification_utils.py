#!/usr/bin/env python

import pandas as pd

def write_featurefile(sm, score=".", group="gene_id", phase="."):
    with open(sm.input[0], 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for line in fhin:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            items = line.split("\t")
            contig = items[0]
            source = items[1]
            method = items[2]
            start = items[3]
            stop = items[4]
            _strand = items[6]
            geneid = items[8]
            if _strand in ['+1', '1', '+']:
                strand = '+'
            else:
                strand = '-'
            gene_id = "{} {}\n".format(group, geneid.split(";")[0].split("=")[-1])
            fhout.write("\t".join([contig, source, method, start, stop, score,
                                   strand, phase, gene_id]))
    return


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


def get_readlength(f):


def normalize_featurecount(sm):
    df = pd.read_csv(sm.input[0], skiprows=1, sep="\t")
    sample_unit = sm.params.s
    df.columns = list(df.columns)[0:-1] + [sample_unit]
    # Get average mapped read length
    readlength = get_readlength(sm.input[1])

    # Perform normalization
    df = calculate_tpm(df, readlength)

    df_tpm = df.iloc[:, [0, -1]]
    df_tpm.columns = ["gene_id", sample_unit]
    df_tpm.to_csv(sm.output[0], sep="\t", index=False)

    df_raw = df.iloc[:, [0, 6]]
    df_raw.columns = ["gene_id", sample_unit]
    df_raw.to_csv(sm.output[1], sep="\t", index=False)


def main(sm):
    toolbox = {"write_featurefile": write_featurefile,
               "normalize_featurecount": normalize_featurecount}

    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
