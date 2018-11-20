#!/usr/bin/env python

import sys, os, csv, argparse, logging
from KrakenResults import *

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filteredFile', dest='filteredFile', help="the file with the filtered classifier results",
                        type=str, default=None, required=True)
    parser.add_argument('--readAssignmentFile', dest='readAssignmentFile', help="the file with the kraken results",
                        type=str, default=None, required=True)
    parser.add_argument('--outputFile', dest='outputFile', help="the file with the addition of statistics", type=str,
                        default=None, required=True)
    parser.add_argument('--outputDir', dest='outputDir', help="the folder in which the reads are written", type=str,
                        default=None, required=True)
    parser.add_argument('--ncbiTaxonomyFolder', dest='ncbiTaxonomyFolder',
                        help="Where to store the ncbi taxonomy database file (ete3 format)", type=str, required=True)
    parser.add_argument('--printStat', dest='printStat', action='store_false',
                        help="should the statistics of read length and count be printed?", default=True)
    parser.add_argument('--classifier', dest='classifier', default='kraken',
                        help="[kraken/centrifuge] Type of classifier used. Default is kraken")
    parser.add_argument('-v', dest='verbose', action='store_true', default=False)
    args = parser.parse_args()
    outputReadsTo = args.outputDir
    # for each of the potential pathogens found in the previous filtering step
    # calculcate some statistics and pull out the reads assigned to this taxa

    # read in the NCBI taxonomy tree
    logging.info("Getting taxonomy database in {}".format(args.ncbiTaxonomyFolder))
    NCBITaxTree = make_ncbi_taxa(args.ncbiTaxonomyFolder)

    # read in the results from the pathogen, species, read depth and coverage filtering step
    potentials = {}
    genusInsteadOfSpecies = {}
    with open(args.filteredFile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i == 0:
                continue
            taxID = row[0]
            potentials[taxID] = row # remember the whole row of 'annotations'
            if row[6] != "S":
                genusInsteadOfSpecies[taxID] = True

    # read in the transient results file (eg the one that for each read specifies which tax ids were 'found')
    logging.info("Reading {} results from {}".format(args.classifier, args.readAssignmentFile))
    if args.classifier == "kraken":
        readAssignments = KrakenReadLevelClassifications(args.readAssignmentFile,
                args.printStat, potentials, NCBITaxTree)
    elif args.classifier == "centrifuge":
        readAssignments = CentrifugeReadLevelClassifications(args.readAssignmentFile,
                args.printStat, potentials, NCBITaxTree)
    else:
        sys.exit("UNSUPPORTED CLASSIFIER {}\n".format(args.classifier))

    # for each taxa make a file with the read assignments
    if not os.path.isdir(outputReadsTo):
        os.makedirs(outputReadsTo)
    if not os.path.isdir(outputReadsTo+"/species"):
        os.makedirs(outputReadsTo+"/species")
    if not os.path.isdir(outputReadsTo+"/genus"):
        os.makedirs(outputReadsTo+"/genus")
    with open(args.outputFile, 'w') as outFile:
        header ="NCBITaxID\tReadDepthForThisTax\tReadDepthForThisOrLower\tTaxonomy Scientific Name\tGenomeSize\tPossibleCoverage\tRankCode"
        if args.printStat:
            header+="\t"+"MedianSpeciesPerRead\t25Percentile\t75Percentile\tMedianReadLengthForSpecies"

        header+="\n"
        outFile.write(header)
        #Loop over the taxids for which reads have been assigned. Note that this is not
        #necessarily the same as in the krakens sample report, since the latter also
        #contains the taxons all the way up to root
        for taxID in [x for x in potentials if x in readAssignments.species.keys()]:
            # loop through all reads assigned to this species and figure out the median (and 25 and 75 percentile) of how many species were assigned each time
            out = "\t".join(potentials[taxID])
            if args.printStat:
                s = readAssignments.getStatistics(taxID)
                out+="\t{}\t{}\t{}\t{}".format(s[1],s[0],s[2],s[3])
            outFile.write(out+"\n")

            # print all the read identifiers for this species to file
            folder = "species"
            if(taxID in genusInsteadOfSpecies):
                folder = "genus"
            with open(outputReadsTo+"/"+folder+"/"+taxID+".readList", 'w') as speciesFile:
                for r in readAssignments.species[taxID]:
                    speciesFile.write(r+"\n")
                # for a species, I should also go one step up in the NCBI TaxTree and pull out all reads that map to that level...
                if taxID in genusInsteadOfSpecies:
                    lineage = NCBITaxTree.get_lineage(taxID)
                    try:
                        parent = lineage[-2]
                    except IndexError:
                        continue
                    if parent in readAssignments.species:
                        for r in readAssignments.species[parent]:
                            speciesFile.write(r+"\n")


if __name__ == "__main__":
    main()
