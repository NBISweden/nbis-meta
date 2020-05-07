#!/usr/bin/env python

import pandas as pd
df = pd.read_csv(snakemake.input, sep="\t")
dfc = df.groupby(["tRNA_type","Bin_Id"]).count().reset_index().loc[:,["tRNA_type","tRNA#","Bin_Id"]]
table = dfc.pivot_table(columns="tRNA_type", index="Bin_Id")["tRNA#"].fillna(0)
table.index.name = table.columns.name = ""
table.to_csv(snakemake.output[0], sep="\t", index=True)

total = {}
for m in table.index:
    c = len(table.loc[m, table.loc[m]>0])
    total[m] = c
table = pd.DataFrame(total, index=["tRNAs"]).T
table.index.name = "Bin_Id"
table.to_csv(snakemake.output[1], sep="\t", index=True)