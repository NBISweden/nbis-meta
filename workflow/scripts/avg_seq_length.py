#!/usr/bin/env python

import pandas as pd

df = pd.read_csv(snakemake.input[0], header=0, sep="\t", index_col=0)
df = df.loc[:, "FastQC_mqc-generalstats-fastqc-avg_sequence_length"]
df.columns = ["read_length"]
df.to_csv(snakemake.output[0], sep="\t", index=True)
