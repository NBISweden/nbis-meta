import sys
from argparse import ArgumentParser

def parse_lines(fh, args):
    for line in fh:
        if line[0] == '#':
            continue
        contig, source, method, start, stop, _, _strand, _, geneid = line.split("\t")
        if args.source:
            source = args.source
        if args.method:
            method = args.method
        if _strand in ['+1', '1', '+']:
            strand = '+'
        else:
            strand = '-'
        gene_id = args.group + " "+ geneid.split(";")[0].split("=")[-1]+"\n"
        sys.stdout.write("\t".join([contig, source, method, start, stop, args.score, strand, args.phase, gene_id]))

def main():
    parser = ArgumentParser('''Parser script to produce featurecount compatible GFF2 file from prodigal orf call''')
    parser.add_argument("gff", nargs=1)
    parser.add_argument("-s", "--source")
    parser.add_argument("-m", "--method")
    parser.add_argument("--score", default=".")
    parser.add_argument("-p", "--phase", default=".")
    parser.add_argument("-g", "--group", default="gene_id")

    args = parser.parse_args()

    with open(args.gff[0], 'r') as fh:
        parse_lines(fh, args)

if __name__ == '__main__':
    main()