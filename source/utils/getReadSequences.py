#!/usr/bin/python

import argparse, os, glob, sys

parser = argparse.ArgumentParser()
parser.add_argument('--folder', dest='readListfolder', help="the folder containing one file per potential pathogen of read identifiers", type=str, default=None, required=True);
parser.add_argument('--fastQFile', dest='fastQFile', help="the .fastQ file of reads", type=str, default=None, required=True);
parser.add_argument('--fakeOut', dest='fakeOut', help="do not use", type=str, default=None, required=True);
args = parser.parse_args()

fileOfAllReadIdentifiers = args.readListfolder + "/_allMappedReadIdentifiers"
fileOfAllReadSequences = fileOfAllReadIdentifiers + ".fastq"
cmd = "cat "+args.readListfolder + "/* | sort -u > "+fileOfAllReadIdentifiers
os.system(cmd) # take all unique read identifiers
cmd = "seqtk subseq "+args.fastQFile+" "+fileOfAllReadIdentifiers+" > "+fileOfAllReadSequences
os.system(cmd) # pull out the sequences for all read identifiers

# now for each potential pathogen pull out only those read sequences
files = glob.glob(args.readListfolder+"/*.readList")
for f in files:
    cmd = "seqtk subseq "+fileOfAllReadSequences+ " "+f+" > "+f + ".fastq"
    os.system(cmd)
