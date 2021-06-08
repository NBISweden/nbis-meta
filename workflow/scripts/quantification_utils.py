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


def clean_featurecount(sm):
    """
    This cleans the featureCounts output table from format:
    # Program:featureCounts v2.0.0; Command:"featureCounts" "-a" "etc."
    Geneid  Chr         Start   End     Strand  Length  path/to/bam/sample.bam
    1_1     k141_7581   1       459     -       459     1
    2_1     k141_0      469     714     +       246     2

    To format:
    gene_id         Length  sample
    k141_7581_1     459     1
    k141_0_1        246     1
    """
    df = pd.read_csv(sm.input[0], comment="#", sep="\t")
    # Extract gene number and combine with contig id
    df["gene_num"] = [x[1] for x in df.Geneid.str.split("_")]
    df.set_index(df.Chr.map(str) + "_" + df.gene_num, inplace=True)
    df.drop("gene_num", axis=1, inplace=True)
    df.index.name = 'gene_id'
    # Set sample and unit name from wildcards
    sample_unit = "{sample}_{unit}".format(sample=sm.wildcards.sample,
                                           unit=sm.wildcards.unit)
    df.columns = list(df.columns)[0:-1] + [sample_unit]
    # Extract length and counts
    df = df.loc[:, ["Length", sample_unit]]
    df.to_csv(sm.output[0], sep="\t")


def aggregate_featurecount(sm):
    """
    Aggregates cleaned featureCounts tables into one table per assembly

    :param sm: snakemake object
    :return:
    """
    df = pd.DataFrame()
    lmap = {}
    for f in sm.input:
        _df = pd.read_csv(f, sep="\t", index_col=0)
        lmap.update(_df.to_dict()["Length"])
        _df.drop("Length", axis=1, inplace=True)
        df = pd.merge(df, _df, right_index=True, left_index=True, how="outer")
    counts = pd.merge(df, pd.DataFrame(lmap, index=["Length"]).T,
                      left_index=True, right_index=True)
    counts.to_csv(sm.output[0], sep="\t")


def process_and_sum(q_df, annot_df):
    # Merge annotations and abundance
    # keep ORFs without annotation as "Unclassified"
    annot_q_df = pd.merge(annot_df, q_df, left_index=True, right_index=True,
                          how="right")
    annot_q_df.fillna("Unclassified", inplace=True)
    feature_cols = annot_df.columns
    annot_q_sum = annot_q_df.groupby(list(feature_cols)).sum().reset_index()
    annot_q_sum.set_index(feature_cols[0], inplace=True)
    return annot_q_sum


def sum_to_features(abundance, parsed):
    parsed_df = pd.read_csv(parsed, index_col=0, sep="\t")
    abundance_df = pd.read_csv(abundance, index_col=0, sep="\t")
    abundance_df.drop("Length", axis=1, inplace=True, errors="ignore")
    feature_sum = process_and_sum(abundance_df, parsed_df)
    return feature_sum


def count_features(sm):
    """
    Counts reads mapped to features such as KOs, PFAMs etc.

    :param sm:
    :return:
    """
    feature_sum = sum_to_features(sm.input.abund, sm.input.annot)
    feature_sum.to_csv(sm.output[0], sep="\t")

def sum_to_taxa(sm):
    """
    Takes taxonomic assignments per orf and abundance values (tpm or raw) and
    returns abundance values summed to unique combinations of ranks

    :param sm: snakemake object
    :return:
    """
    header = ["protein", "superkingdom", "phylum", "class", "order", "family",
              "genus", "species"]
    df = pd.read_csv(sm.input.tax[0], sep="\t", index_col=0, header=None,
                     names=header)
    abund_df = pd.read_csv(sm.input.abund, header=0, index_col=0, sep="\t")
    # Remove length column
    abund_df.drop("Length", axis=1, inplace=True, errors="ignore")
    taxa_abund = pd.merge(df, abund_df, right_index=True, left_index=True)
    taxa_abund_sum = taxa_abund.groupby(header[1:]).sum().reset_index()
    taxa_abund_sum.to_csv(sm.output[0], sep="\t", index=False)


def main(sm):
    toolbox = {"write_featurefile": write_featurefile,
               "clean_featurecount": clean_featurecount,
               "aggregate_featurecount": aggregate_featurecount,
               "count_features": count_features,
               "sum_to_taxa": sum_to_taxa}

    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
