#!/usr/bin/env python

import pandas as pd


def parse_pfam(sm):
    annot = pd.read_csv(sm.input[0], comment="#", header=None, sep=" +",
                        usecols=[0, 5, 7, 14], engine="python",
                        names=["orf", "pfam", "pfam_type", "pfam_clan"])
    clans = pd.read_csv(sm.input[1], header=None, names=["clan", "clan_name"],
                        usecols=[0, 3], sep="\t")
    info = pd.read_csv(sm.input[2], header=None,
                       names=["pfam", "clan", "pfam_name"], usecols=[0, 1, 4],
                       sep="\t")
    # Strip suffix for pfams
    annot.loc[:, "pfam"] = [x.split(".")[0] for x in annot.pfam]
    # Select unique orf->pfam mappings
    # TODO: This masks multiple occurrences of domains on the same orf. Figure out if this is wanted or not.
    # Merge with pfam info and clan info
    annot = annot.groupby(["orf", "pfam"]).first().reset_index()
    annot = pd.merge(annot, info, left_on="pfam", right_on="pfam")
    annot = pd.merge(annot, clans, left_on="clan", right_on="clan", how="left")
    annot.fillna("No_clan", inplace=True)
    annot = annot.loc[:, ["orf", "pfam", "pfam_name", "clan", "clan_name"]]
    annot.sort_values("orf", inplace=True)
    # Write to file
    annot.to_csv(sm.output[0], sep="\t", index=False)


def main(sm):
    toolbox = {"parse_pfam": parse_pfam}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
