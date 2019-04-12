import sys, os, csv, re, argparse, operator

parser = argparse.ArgumentParser()
parser.add_argument('--data', dest='inFile', type=str, default=None, required=True);
parser.add_argument('--alignmentFolder', dest='alignmentFolder', type=str, default=None, required=True);
parser.add_argument('--output', dest='outputFile', type=str, default=None, required=True);
args = parser.parse_args()


def main():
    # loop through the file and add a quality assessment for each hit...
    df = pd.read_csv(args.inFile, sep="\t")    # "NCBITaxID\tTaxonomy\tIndividual\tRun\tCoverage\tNrReads\tAvgReadLength\tNormalisedReadCount\tAncientStatus\tAssemblyLevel\tNumberOfEntriesInFasta\tSequenceDepth\n"
    with open(args.outputFile, 'w') as fo:
        fo.write("NCBITaxID\tTaxonomy\tIndividual\tRun\tScore\n")
        for index, row in df.iterrows():
            score = 0
            if row['Coverage'] > 0.8:
                score = score + 1
            elif row['Coverage'] > 0.5:
                score = score + 0.6
            elif row['Coverage'] > 0.3:
                score = score + 0.1
            if row['NormalisedReadCount'] > 1000:
                score = score + 1
            elif row['NormalisedReadCount'] > 500:
                score = score + 0.6
            elif row['NormalisedReadCount'] > 100:
                score = score + 0.1
            if row['AssemblyLevel'] == "Complete Genome":
                score = score + 1
            elif row['AssemblyLevel'] == "Chromosome":
                score = score + 0.6
            elif row['AssemblyLevel'] == "Scaffold":
                score = score + 0.2
            elif row['AssemblyLevel'] == "Contig":
                score = score + 0.1
            if row['NumberOfEntriesInFasta'] < 2:
                score = score + 1
            elif row['NumberOfEntriesInFasta'] < 4:
                score = score + 0.6
            elif row['NumberOfEntriesInFasta'] < 5:
                score = score + 0.5
            if row['SequenceDepth'] > 1000:
                score = score + 1
            elif row['SequenceDepth'] > 500:
                score = score + 0.5
            bam=os.path.join(args.alignmentFolder,row['Individual'],row['Run'],row['NCBITaxID']+".bam")
            score = score + fullAlignment(bam)
            fo.write(row['NCBITaxID']+"\t"+row['Taxonomy']+"\t"+row['Individual']+"\t"+row['Run']+"\t"+str(score))


# for a given hit, check how many of the reads remaining after I align to the full database...
def fullAlignment(file):
    return(1)



if __name__ == "__main__":
    main()
