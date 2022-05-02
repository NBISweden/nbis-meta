#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import sys
from subprocess import run
import math


def count_seqs(f):
    """
    Uses subprocess to call 'grep' and count number of entries in fasta
    file

    :param f: Input fasta file
    :return: integer of sequences in file
    """
    p = run(['grep', '-c', '>', f],
            capture_output=True)
    return int(p.stdout.decode().rstrip())


def write_files(f, n_files, n_seqs, prefix="split", outdir="."):
    i = 1
    fhout = open(f"{outdir}/{prefix}_{i}.fasta", 'w')
    sys.stderr.write(f"Writing to files under {outdir}:\n")
    with open(f, 'r') as fhin:
        for j, record in enumerate(parse(f, "fasta"), start=1):
            fhout.write(f">{record.description}\n{record.seq}\n")
            if j == n_seqs:
                break
            if j % n_files == 0:
                i+=1
                fhout.close()
                fhout = open(f"{outdir}/{prefix}_{i}-of-{n_files}.fasta", 'w')
    fhout.close()


def main(args):
    if not args.n_seqs:
        n_seqs = count_seqs(args.infile)
    else:
        n_seqs = args.n_seqs
    seqs_per_file = math.ceil(n_seqs/args.n_files)
    sys.stderr.write(f"Files: {args.n_files}\n")
    sys.stderr.write(f"Sequences: {n_seqs}\n")
    sys.stderr.write(f"Sequences_per_file: {seqs_per_file}\n")
    write_files(args.infile, seqs_per_file, n_seqs, prefix=args.prefix, outdir=args.outdir)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Input file")
    parser.add_argument("-n", "--n_files", type=int, default=10,
                        help="Number of files to split into")
    parser.add_argument("-N", "--n_seqs", type=int,
                        help="Number of sequences in input file")
    parser.add_argument("-o", "--outdir", type=str, default=".",
                        help="Output directory for split files")
    parser.add_argument("--prefix", type=str, default="split",
                        help="Prefix for split files")
    args = parser.parse_args()
    main(args)
