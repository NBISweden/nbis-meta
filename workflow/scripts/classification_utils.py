#!/usr/bin/env python

import pandas as pd


def metaphlan2krona(sm):
    """
    Creates a krona-input table from metaphlan (standard format) results
    :param sm: snakemake object
    :return:
    """
    df = pd.read_csv(sm.input[0], sep="\t", header=None, comment="#",
                     index_col=0, usecols=[0, 1, 2],
                     names=["lineage", "taxids", "%"])
    # Extract species level
    df_sp = df.loc[df.index.str.contains("s__")]
    # Create new dataframe with taxids \t %
    data = dict(zip([x.split("|")[-1] for x in df_sp.taxids], df_sp["%"]))
    df_out = pd.DataFrame(data, index=["%"]).T
    df_out.to_csv(sm.output[0], sep="\t", index=True, header=False)


def main(sm):
    toolbox = {"metaphlan2krona_table": metaphlan2krona}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
