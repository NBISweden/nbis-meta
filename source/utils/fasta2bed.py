from Bio.SeqIO import parse
from argparse import ArgumentParser
import sys

def write_regions(f,fh, saf=False):
    for record in parse(f, "fasta"):
        if saf:
            fh.write("{}\t{}\t{}\t{}\t{}\n".format(record.id,record.id,0,len(record),"+"))
        else:
            fh.write("{}\t{}\t{}\n".format(record.id,0,len(record)))

def regions(args):
    if args.outfile:
        with open(args.outfile, 'w') as fh:
            write_regions(args.infile, fh, args.saf)
    else:
        write_regions(args.infile, sys.stdout, args.saf)

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
                        help="Input fasta file")
    parser.add_argument("-o", "--outfile", type=str,
                        help="Output bed file")
    parser.add_argument("--saf", action="store_true",
                        help="If specified, write SAF (simplified annotation format) output")
    args = parser.parse_args()

    regions(args)

if __name__ == '__main__':
    main()
