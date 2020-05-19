#!/usr/bin/env python

import pandas as pd


def write_featurefile(sm, score=".", group="gene_id", phase="."):
    """
    Takes a prodigal GFF file and turns it into a file for use with
    featureCounts
    :param sm: snakemake object
    :param score: dummy score column
    :param group: identifier for featurecounts
    :param phase: dummy phase column
    :return:
    """
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


def calculate_tpm(df, readlength, fhlog):
    """
    Calculates transcripts per million normalized values for a sample based
    on featureCounts output
    :param df: pandas DataFrame as read from featureCounts output
    :param readlength: average length of mapped reads
    :param fhlog: loghandle
    :return: pandas DataFrame with normalized values per gene
    """
    sampleName = df.columns[-1]
    # 1. Calculate t for sample
    # t = (reads_mapped_to_gene * read_length) / length_of_gene
    # Multiply gene counts with read length,
    # then divide by the Length column (in kb)
    fhlog.write("Normalizing by read length and gene length\n")
    t = df[sampleName].multiply(readlength).div(df["Length"].div(1000))
    df = df.assign(t=pd.Series(t, index=df.index))
    # 2. Calculate T
    # T = sum(t)
    fhlog.write("Calculating sum of normalized values\n")
    T = df["t"].sum()
    # 3. Calculate TPM
    # TPM = t*10^6 / T
    fhlog.write("Normalization factor T is {}\n".format(T))
    fhlog.write("Calculating TPM\n")
    TPM = (df["t"].multiply(1000000)).div(T)
    df = df.assign(TPM=pd.Series(TPM, index=df.index))
    return df


def get_readlength(f):
    """
    Reads samtools flagstat file and extract average mapped read length
    :param f: samtools flagstat file
    :return:
    """
    with open(f, 'r') as fhin:
        for line in fhin:
            line = line.rstrip()
            if line.startswith("SN") and "average length:" in line:
                length = line.rsplit()[-1]
    return float(length)


def normalize_featurecount(sm):
    """
    Master rule for running normalization
    :param sm: snakemake object
    :return:
    """
    df = pd.read_csv(sm.input[0], skiprows=1, sep="\t")
    sample_unit = "{sample}_{unit}".format(sample=sm.wildcards.sample,
                                           unit=sm.wildcards.unit)
    df.columns = list(df.columns)[0:-1] + [sample_unit]
    # Get average mapped read length
    readlength = get_readlength(sm.input[1])

    # Perform normalization
    with open(sm.log[0], 'w') as fhlog:
        df = calculate_tpm(df, readlength, fhlog)

    df_tpm = df.iloc[:, [0, -1]]
    df_tpm.columns = ["gene_id", sample_unit]
    df_tpm.to_csv(sm.output[0], sep="\t", index=False)

    df_raw = df.iloc[:, [0, 6]]
    df_raw.columns = ["gene_id", sample_unit]
    df_raw.to_csv(sm.output[1], sep="\t", index=False)


def merge_files(files, gff_df):
    """
    Merges abundance tables from several samples
    :param files: list of files
    :param gff_df: pandas DataFrame with orf ids in the format <contig>_<orfnum>
    :return: merged pandas DataFrame
    """
    df = pd.DataFrame()
    for f in files:
        _df = pd.read_csv(f, index_col=0, sep="\t")
        df = pd.concat([df, _df], axis=1)
    df = pd.merge(df, gff_df, left_index=True, right_on="gene_id")
    df.drop("gene_id", axis=1, inplace=True)
    df.set_index("orf", inplace=True)
    return df


def aggregate_featurecount(sm):
    """
    Aggregates normalized and raw tables per sample into one table per assembly

    :param sm: snakemake object
    :return:
    """
    gff_df = pd.read_csv(sm.input.gff_file, header=None, usecols=[0, 8],
                         names=["contig", "gene"], sep="\t")
    gff_df = gff_df.assign(
        gene_id=pd.Series([x.replace("gene_id ", "") for x in gff_df.gene],
                          index=gff_df.index))
    gff_df = gff_df.assign(
        suffix=pd.Series([x.split(" ")[-1].split("_")[-1] for x in gff_df.gene],
                         index=gff_df.index))
    gff_df = gff_df.assign(
        orf=pd.Series(gff_df.contig + "_" + gff_df.suffix, index=gff_df.index))
    gff_df = gff_df[["orf", "gene_id"]]

    raw_df = merge_files(sm.input.raw_files, gff_df)
    tpm_df = merge_files(sm.input.tpm_files, gff_df)
    raw_df.to_csv(sm.output.raw, sep="\t")
    tpm_df.to_csv(sm.output.tpm, sep="\t")


def sum_to_taxa(sm):
    """
    Takes taxonomic assignments per orf and abundance values (tpm or raw) and
    returns abundance values summed to unique combinations of ranks

    :param sm: snakemake object
    :return:
    """
    header = ["protein", "superkingdom", "phylum", "class", "order", "family",
              "genus", "species"]
    df = pd.read_csv(sm.input.tax, sep="\t", index_col=0, header=None,
                     names=header)
    abund_df = pd.read_csv(sm.input.abund, header=0, index_col=0, sep="\t")
    taxa_abund = pd.merge(df, abund_df, right_index=True, left_index=True)
    taxa_abund_sum = taxa_abund.groupby(header[1:]).sum().reset_index()
    taxa_abund_sum.to_csv(sm.output[0], sep="\t", index=False)


def sum_to_rgi(sm):
    """
    Takes Resistance gene identifier output and abundance info per orf and sums
    values to gene family level

    :param sm: snakemake object
    :return:
    """
    annot = pd.read_csv(sm.input.annot, sep="\t", header=0, index_col=0,
                        usecols=[0, 16])
    # Rename index for annotations to remove text after whitespace
    annot.rename(index=lambda x: x.split(" ")[0], inplace=True)
    abund = pd.read_csv(sm.input.abund, sep="\t", header=0, index_col=0)
    df = pd.merge(annot, abund, left_index=True, right_index=True)
    # Sum to Gene family
    dfsum = df.groupby("AMR Gene Family").sum()
    dfsum.to_csv(sm.output[0], sep="\t", index=True, header=True)

def main(sm):
    toolbox = {"write_featurefile": write_featurefile,
               "normalize_featurecount": normalize_featurecount,
               "aggregate_featurecount": aggregate_featurecount,
               "sum_to_taxa": sum_to_taxa,
               "sum_to_rgi": sum_to_rgi}

    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
