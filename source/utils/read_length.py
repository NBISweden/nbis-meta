#!/usr/bin/env python
import argparse
import subprocess


def main(args):
    assert len(args.input_fastqs) % 2 == 0
    if len(args.sample_names)*2 != len(args.input_fastqs):
        print(args.sample_names, args.input_fastqs)
    assert len(args.sample_names)*2 == len(args.input_fastqs)

    sample_names_and_fastqs = zip(args.sample_names, zip(args.input_fastqs[0::2], args.input_fastqs[1::2]))
    for sample_name, fastq_tuple in sample_names_and_fastqs:
        R1, R2 = fastq_tuple

        # Use a random 1 % of the reads
        proc = subprocess.Popen(['seqtk','sample', R1, '0.01'], stdout=subprocess.PIPE)
        R1_total_length = 0
        for i, line in enumerate(proc.stdout):
            if i % 4 == 1:
                line.strip()
                R1_total_length+=len(line)
        R1_nr_reads = (i+1)/4

        R2_total_length = 0
        # Use a random 1 % of the reads
        proc = subprocess.Popen(['seqtk','sample', R2, '0.01'], stdout=subprocess.PIPE)
        for i, line in enumerate(proc.stdout):
            if i % 4 == 1:
                line.strip()
                R2_total_length+=len(line)

        R2_nr_reads = (i+1)/4
        avg_read_length = (R1_total_length + R2_total_length) / float(R1_nr_reads + R2_nr_reads)

        print("{0}\t{1:0.4f}".format(sample_name, avg_read_length))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fastqs", nargs='*')
    parser.add_argument("--sample_names", nargs='*')

    args = parser.parse_args()
    main(args)
