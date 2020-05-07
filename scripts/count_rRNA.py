#!/usr/bin/env python

import pandas as pd

df = pd.read_csv(snakemake.input, sep="\t", usecols=[0, 2, 8], header=None,
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
table.to_csv(snakemake.output, sep="\t", index=True)
