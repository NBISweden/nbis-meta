#!/usr/bin/env python

from argparse import ArgumentParser
from ete3 import NCBITaxa
import sys
import pandas as pd

def read_centrifuge_report(reportfile):
    df = pd.read_table(reportfile, header=0)
    tmp = df.copy(deep=True)
    genomesizes = tmp.to_dict()["genomeSize"]
    return df, genomesizes


def read_centrifuge_result(args):
    taxa_counts = {}
    with open(args.infile, 'r') as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            line = line.rstrip()
            readID, seqID, taxID, score, secondBestScore, hitLength, queryLength, numMatches = line.split("\t")
            if args.unique and int(numMatches) > 1:
                continue

            if args.min_length and int(hitLength) < args.min_length:
                continue

            if args.min_score and int(score) < args.min_score:
                continue

            if not taxID in taxa_counts.keys():
                taxa_counts[taxID] = {"count": 0, "aligned": 0}

            count = 1/float(numMatches)
            aligned = float(hitLength)/float(numMatches)
            taxa_counts[taxID]["count"] += count
            taxa_counts[taxID]["aligned"] += aligned
    return taxa_counts


def generate_ete3db(taxdbfile):
    if not taxdbfile:
        taxdb = NCBITaxa()
    else:
        import os
        import subprocess
        dirname = os.path.dirname(taxdbfile)
        if dirname != "":
            cmd = "mkdir -p "+dirname
            p1 = subprocess.Popen(cmd, shell=True, stdin=None)
            p1.wait()
        cmd = "touch "+taxdbfile
        p2 = subprocess.Popen(cmd, shell=True, stdin=None)
        p2.wait()
        taxdb = NCBITaxa(taxdbfile)
    return taxdb


def get_lineage_names(lineage, taxdb):
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    names = {}
    for rank in ranks:
        names[rank] = ""
    for id, rank in lineage.items():
        if rank not in ranks:
            continue
        name = taxdb.translate_to_names([id])[0]
        names[rank] = name
    name_string = ""
    for rank in ranks:
        name_string+="\t{}".format(names[rank])
    return name_string


def generate_krona_table(taxa_counts, args, genomesizes):
    taxdb = generate_ete3db(args.taxdb)
    lineage_dict = {}
    with open(args.kronaout, 'w') as fhout:
        for taxid in taxa_counts.keys():
            try:
                name_string = lineage_dict[taxid]
            except KeyError:
                try:
                    lineage = taxdb.get_rank(taxdb.get_lineage(taxid))
                    name_string = get_lineage_names(lineage, taxdb)
                except ValueError:
                    name_string = "\t" + "\t".join(["Unknown"]*7)
                lineage_dict[taxid] = name_string
            count = taxa_counts[taxid]["count"]
            aligned = taxa_counts[taxid]["aligned"]
            if args.normalize and genomesizes:
                try:
                    genomesize = genomesizes[int(taxid)]
                except KeyError:
                    continue
                # Normalized is aligned bases per kb of genome
                try:
                    norm = aligned / float(genomesize) * 1000
                except ZeroDivisionError:
                    norm = 0
                fhout.write("{}{}\n".format(norm, name_string))
            else:
                fhout.write("{}{}\n".format(count, name_string))


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
                        help="Centrifuge output.")
    parser.add_argument("-r", "--reportfile", type=str,
                        help="Centrifuge report file. Used to normalize to genome sizes")
    parser.add_argument("--kronaout", type=str,
                        help="Produce output ready to import to Krona using ktImportText")
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize read counts by genome size")
    parser.add_argument("--unique", action="store_true",
                        help="Only count reads mapping uniquely")
    parser.add_argument("--min_score", type=int,
                        help="Require a minimum score for reads to be counted")
    parser.add_argument("--min_length", type=int,
                        help="Require a minimum alignment length to the read")
    parser.add_argument("--taxdb", type=str,
                        help="Ete3 sqlite database file or path to one that will be created.")
    #parser.add_argument("--ranksplit", action="store_true",
    #                    help="")

    args = parser.parse_args()

    sys.stderr.write("Reading centrifuge results\n")
    taxa_counts = read_centrifuge_result(args)
    sys.stderr.write("{} taxa parsed\n".format(len(taxa_counts)))

    if args.reportfile:
        reportdf, genomesizes = read_centrifuge_report(args.reportfile)
    else:
        genomesizes = False

    if args.kronaout:
        generate_krona_table(taxa_counts, args, genomesizes)


if __name__ == '__main__':
    main()