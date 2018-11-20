#!/usr/bin/env python

import pandas as pd
from ete3 import NCBITaxa
from argparse import ArgumentParser
import sys


def make_krona_table(f, db):
    if not db:
        ncbi_taxa = NCBITaxa()
    else:
        ncbi_taxa = NCBITaxa(db)
    krona_table = pd.DataFrame(columns = ["abundance","superkingdom","phylum","class","order","family","genus",
                                          "species","leaf"])
    one_letter_ranks = {"D": "superkingdom", "P": "phylum", "C": "class", "O": "order", "F": "family", "G": "genus",
                        "S": "species"}
    df = pd.read_table(f, header=None, names = ["clade_percent", "clade_reads", "reads", "rank", "taxid", "name"])
    df = df.loc[df.reads > 0]
    for j, i in enumerate(df.index):
        r = df.loc[i]
        taxid = r["taxid"]
        reads = r["reads"]
        name = r["name"]
        one_letter_rank = r["rank"]
        if one_letter_rank == "-":
            rank = ncbi_taxa.get_rank([taxid])[taxid]
            try:
                parent_taxid = ncbi_taxa.get_lineage(taxid)[-2]
            except IndexError:
                parent_taxid = taxid
            parent_rank = ncbi_taxa.get_rank([parent_taxid])[parent_taxid]
            if rank == "no rank" and parent_rank == "species":
                rank = "leaf"
            else:
                continue
        elif one_letter_rank == "U":
            rank = "unclassified"
        else:
            try:
                rank = one_letter_ranks[one_letter_rank]
                #TODO: Shouldn't be too many reads mapped directly to ranks not in the krona table, but check eventually
            except KeyError:
                continue
        res = {"abundance": reads, "superkingdom": "", "phylum": "", "class": "", "order": "", "family": "",
               "genus": "", "species": "", "leaf": ""}
        if rank != "unclassified":
            rank_dict = ncbi_taxa.get_rank(ncbi_taxa.get_lineage(taxid))
            name_dict = ncbi_taxa.get_taxid_translator(ncbi_taxa.get_lineage(taxid))
            for dict_taxid, dict_rank in rank_dict.items():
                if dict_rank in res.keys():
                    rank_name = name_dict[dict_taxid]
                    res[dict_rank] = rank_name
            if not rank in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
                res["leaf"] = name
        _df = pd.DataFrame(res, index=[j])[krona_table.columns]
        krona_table = pd.concat([krona_table, _df])
    return krona_table


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
                        help="Centrifuge output file")
    parser.add_argument("-t", "--taxdb", default=False,
                        help="Specific ete3 NCBI sqlite database. Leave empty to use (or initialize) the one in your home directory.")
    args = parser.parse_args()

    krona_table = make_krona_table(args.input, args.taxdb)
    krona_table.to_csv(sys.stdout, sep="\t", index=False, header=False)

if __name__ == '__main__':
    main()