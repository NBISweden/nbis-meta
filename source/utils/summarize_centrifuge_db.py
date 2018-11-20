#!/usr/bin/env python

from argparse import ArgumentParser
from ete3 import NCBITaxa
import sys

def make_seq2sizemap(f):
    seq2sizemap = {}
    sys.stderr.write("Reading sequence information from {}\n".format(f))
    total_size = 0
    with open(f, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            try:
                _, identifier, size = line.split("\t")
            except ValueError:
                continue
            seq = identifier.split(" ")[0]
            size = int(size)
            seq2sizemap[seq] = size
            total_size+=size
    sys.stderr.write("Read {} sequences. Total size = {}bp\n".format(len(seq2sizemap), total_size))
    return seq2sizemap


def summarize_to_ranks(seq2sizemap, mapfile, summary_rank, update):
    ncbitaxa = NCBITaxa()
    if update:
        ncbitaxa.update_taxonomy_database()
    summary = {}
    sys.stderr.write("Reading seq2taxid mapfile and summarizing at rank {}\n".format(summary_rank))
    with open(mapfile, 'r') as fh:
        for i, line in enumerate(fh, start=1):
            if i%1000 == 0:
                sys.stderr.write("Read {} lines...\n".format(i))
            line = line.rstrip()
            try:
                seqid, taxid = line.split("\t")
            except ValueError:
                seqid, taxid = line.split()
            taxid = int(taxid)
            try:
                size = seq2sizemap[seqid]
            except KeyError:
                sys.stderr.write("WARNING: Sequence {} missing from summary\n".format(seqid))
                continue
            try:
                ranks = ncbitaxa.get_rank(ncbitaxa.get_lineage(taxid))
            except ValueError:
                sys.stderr.write("WARNING: Taxid {} missing from db\n".format(taxid))
                continue
            if summary_rank not in ranks.values():
                rank_name = "Unknown"
            else:
                for rank_taxid, rank in ranks.items():
                    if rank == summary_rank:
                        rank_name = ncbitaxa.get_taxid_translator([rank_taxid])[rank_taxid]
            if not rank_name in summary:
                summary[rank_name] = {"size": 0, "seqs": 0}
            summary[rank_name]["size"]+=size
            summary[rank_name]["seqs"]+=1
    return summary


def main():
    parser = ArgumentParser()
    parser.add_argument("--seqid2taxidmap", required=True,
                        help="Sequence id to taxid mapping file.")
    parser.add_argument("--summary", required=True,
                        help="Summary output produced using 'centrifuge-inspect -s'")
    parser.add_argument("--rank", default="superkingdom",
                        help="Summarize sequences and sizes at this rank (default = superkingdom)")
    parser.add_argument("--update", action="store_true",
                        help="Update NCBI taxonomy db before summarizing (default = False)")
    args = parser.parse_args()

    seq2sizemap = make_seq2sizemap(args.summary)

    summary = summarize_to_ranks(seq2sizemap, args.seqid2taxidmap, args.rank, args.update)

    print("{}\t{}\t{}".format("Name","Sequences","TotalSize"))
    for rank_name in sorted(summary.keys()):
        print("{}\t{}\t{}".format(rank_name, summary[rank_name]["seqs"], summary[rank_name]["size"]))

if __name__ == '__main__':
    main()