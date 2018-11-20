#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to calculate TPM values for contigs or genes based on count files

TPM values are defined as in Wagner et al (Theory in Biosciences) 2012. 

      rg x rl x 10^6
TPM = --------------
        flg x T

rg: reads mapped to gene g
rl: read length
flg: feature length (in kb)
T: sum of rgxrl/flg for all genes

Cloned from https://github.com/EnvGen/toolbox/
"""
import pandas as pd, sys, argparse, logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def gene_lengths_from_gff(gff_file, feature_type="CDS", saf=False):
    gene_lengths = {}
    gene_ids = []
    with open(gff_file) as fh:
        for line in fh:
            line = line.rstrip()
            if saf:
                gene_id = line.split("\t")[0]
                gene_lengths[gene_id] = int(line.split("\t")[3]) - int(line.split("\t")[2]) + 1
            else:
                ft = line.split("\t")[2]
                if ft != feature_type:
                    continue
                gene_id = line.split("\t")[-1].replace("gene_id ","")
                gene_lengths[gene_id] = int(line.split("\t")[4])-int(line.split("\t")[3])+1
            gene_ids.append(gene_id)
    return pd.Series(gene_lengths),gene_ids


def main(args):
    logging.info("Reading sample info")
    sample_info = pd.read_table(args.sample_info, header=0, index_col=0, names=['avg_read_len'])
    logging.info("Reading gene lengths from gff")

    gene_lengths,gene_ids = gene_lengths_from_gff(args.gff, args.saf)

    df = pd.DataFrame()
    first = True
    for fn in args.coverage_files:
        # Read counts per gene for sample
        rg = pd.read_table(fn, index_col=0, header=0, compression=args.input_compression)
        sample_name = list(rg.columns)[-1]
        logging.info("Calculating TPM for " + sample_name)
        # Intersect with genes in the gene length file
        rg = rg.loc[list(set(gene_lengths.index).intersection(set(rg.index)))]
        gene_lengths = gene_lengths.loc[list(rg.index)]

        # Average read length for sample
        rl = sample_info.ix[sample_name,'avg_read_len']
        ## Calculate T for sample
        T = rl * rg[sample_name].divide(gene_lengths).sum()

        # Calculate TPM for sample
        tpm = ((1e6*rl)/float(T))*(rg[sample_name].divide(gene_lengths))

        # Create dataframe
        TPM = pd.DataFrame(tpm, columns=[sample_name])

        # Add gene length as the first column
        if first:
            first = False
            df['gene_length'] = gene_lengths
        
        # Concatenate to results
        df = pd.concat([df,TPM],axis=1)

    df.index.name = 'gene_id'
    # Sort output by input
    l = [x for x in gene_ids if x in df.index]
    df = df.loc[l]

    # Write to file
    df.to_csv(sys.stdout, sep='\t')
    logging.info("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('-n', '--sample_names', nargs='*',
    #        help="Sample names, in the same order as coverage_files")
    parser.add_argument('-c', '--coverage_files', nargs='*', 
            help="Coverage files with tab separated values: 'sequence id, count'")
    parser.add_argument('--gff',
            help=("GFF version 2 file"))
    parser.add_argument('--saf', action='store_true',
            help="Input GFF file is actually SAF format")
    parser.add_argument('-i', '--sample_info', 
            help="Tab separated values 'sample_id', 'avg_read_length'")
    parser.add_argument("--input_compression", default=None, choices=[None, 'gzip'], 
            help="Compression type for input coverage files. Default=None, use 'gzip', for gzipped files.")
    args = parser.parse_args()
    main(args)
