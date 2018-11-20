from Bio.SeqIO import parse
import gzip as gz, sys, logging
from argparse import ArgumentParser

def read_list_of_reads(f):
    reads = []
    with open(f, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            reads.append(line)
    return reads


def write_reads(readfile, reads, outfile, verbose):
    if not outfile:
        fh_out = sys.stdout
    else:
        fh_out = open(outfile, 'w')

    if ".gz" in readfile:
        fh_in = gz.open(readfile, 'rt')
    else:
        fh_in = open(readfile, 'r')
    written = 0
    if verbose:
        logging.info("Parsing reads in {}".format(readfile))
    for i, record in enumerate(parse(fh_in, "fastq"),start=1):
        if i%100000 == 0 and verbose:
            logging.info("{} reads parsed, {} reads written...".format(i,written))
        id = record.id
        id_split = id.rsplit("/")[0]
        if id in reads or id_split in reads:
            fh_out.write("{}".format(record.format("fastq")))
            written += 1
    if verbose:
        logging.info("{} reads parsed, {} reads written...Done\n".format(i, written))

def main():
    parser = ArgumentParser()
    parser.add_argument("-r", "--reads", type=str, required = True,
                        help="Reads in (gzipped) fastq format")
    parser.add_argument("-l", "--list", type=str, required=True,
                        help="List of reads to extract")
    parser.add_argument("-o", "--outfile", type=str,
                        help="Write subseqs to outfile (defaults to stdout)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()


    if args.verbose:
        logging.basicConfig(level=logging.INFO,format='%(asctime)s:%(levelname)s:%(message)s')

    logging.info("Getting list of reads from {}".format(args.list))
    reads = read_list_of_reads(args.list)

    write_reads(args.reads, reads, args.outfile, args.verbose)

if __name__ == '__main__':
    main()