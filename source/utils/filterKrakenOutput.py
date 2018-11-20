import sys, os, csv, re, argparse, operator
from KrakenResults import *

parser = argparse.ArgumentParser()
parser.add_argument('--resFile', dest='resFile', help="the file with the classifier results", type=str, default=None, required=True)
parser.add_argument('--whiteListFile', dest='whiteListFile', help="file with taxa of interest", type=str, default=None, required=True)
parser.add_argument('--outputFile', dest='outputFile', help="the file with found taxa", type=str, default=None, required=True)
parser.add_argument('--failedFile', dest='failedFile', help="the file with the kraken results that did not pass all the filters", type=str, default=None, required=True)
parser.add_argument('--genomeSizeFile', dest='genomeSizeFile', help="list of files with genome sizes", type=str, required=True)
parser.add_argument('--readDepthCutoff', dest='readDepthCutoff', help="how many reads have to have been assigned to this species?", type=int, default=10)
parser.add_argument('--readDepthCutoffGenus', dest='readDepthCutoffGenus', help="how many reads have to have been assigned to this genus?", type=int, default=100)
parser.add_argument('--minPotentialCoverage', dest='minPotentialCoverage', help="what should the potential coverage be", type=float, default=0.1)
parser.add_argument('--sortBy', dest='sortBy', help="which field to sort the output on", type=str, default="readDepthThisLevelOnly")
parser.add_argument('-v', dest='verbose', action='store_true', default=False)
args = parser.parse_args()


def main():
    # read in the list of all pathogens, pretty inclusive one at that
    whiteList = WhiteList(args.whiteListFile)

    # read in the genome sizes from NCBI
    genomeSizes = GenomeSizes(args.genomeSizeFile)

    # read in the kraken results file (both the transient one and the constant one) and filter the results
    filterResult(whiteList, genomeSizes, args.resFile)


def filterResult(whiteList, genomeSizes, resFile):

    #if classifier == 'kraken':
    #    # read in the kraken result file - depends on the genomeSizes in order to add genome size in b as well as the potential max coverage
    res = readReportFile(resFile, genomeSizes)
    #elif classifier == 'centrifuge':
        #res = readCentrifugeResultFile(resFile, genomeSizes)

    # show the filtered results
    foundCount = 0
    with open(args.outputFile, 'w') as outFile, open(args.failedFile, 'w') as failedFile:
        outFile.write("NCBITaxID\tReadDepthForThisTax\tReadDepthForThisOrLower\t" \
            "Taxonomy Scientific Name\tGenomeSize\tPossibleCoverage\tRankCode\n")
        failedFile.write("NCBITaxID\tDidNotPassFiltersBecauseOf\n")
        for k in (sorted(res.values(), key=operator.attrgetter(args.sortBy), reverse=True)):
            # see if the current result pass the specified filters, then annotate the result as well - for example with Median Number Of Species Per Read
            filtered = k.passFilters(args.readDepthCutoff, args.readDepthCutoff, args.readDepthCutoffGenus, whiteList, args.minPotentialCoverage)
            if filtered[0]:
                foundCount+=1
                outFile.write(k.printAsTabDelLine())
            else:
                # print the NXBI Tax id along with the reason for failing to the .log file for the step
                failedFile.write(str(k.taxID)+"\t"+filtered[1]+"\n")
            if(k.scientificName == "unclassified"):
                sys.stdout.write(str(k.readDepthThisLevelOnly)+"\tunclassified reads found\n")
            if(k.scientificName == "root"):
                sys.stdout.write(str(k.readDepthThisOrLowerLevel)+"\tclassified reads found with name '"+k.scientificName+"'\n")
    sys.stdout.write("{}\t results were considered\n".format(len(res.keys())))
    sys.stdout.write(str(foundCount)+"\tkraken results passed all the filters\n")

if __name__ == "__main__":
    main()