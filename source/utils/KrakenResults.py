#!/usr/bin/env python

import csv, re, sys, os, logging
import numpy, subprocess

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def make_ncbi_taxa(dbdir=None):
    from ete3 import NCBITaxa
    # Create a ete3 database in your home directory (if no dbfile is specified)
    if not dbdir:
        ncbi_taxa = NCBITaxa()
    # Create a ete3 database in the specified path
    else:
        cmd = "mkdir -p {}".format(dbdir)
        p1 = subprocess.Popen(cmd, shell=True, stdin=None)
        p1.wait()
        dbfile = os.path.join(dbdir,"ncbitax.db")
        cmd = "touch {}".format(dbfile)
        p2 = subprocess.Popen(cmd, shell=True, stdin=None)
        p2.wait()
        ncbi_taxa = NCBITaxa(dbfile)
    return ncbi_taxa


def readReportFile(reportfile, genomeSizes):
    # take a file of kraken results and read in
    results = {}
    with open(reportfile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            taxID = row[4]
            nrReads = int(row[1])
            potentialCoverage = ""
            if genomeSizes is not None:
                genomeSize = genomeSizes.get(taxID)
                potentialCoverage = genomeSizes.coverage(taxID, nrReads)
            else:
                genomeSize = 0
            if potentialCoverage is not "":
                potentialCoverage = float(potentialCoverage)
            kr = KrakenResult(row[0], nrReads, row[2], row[3], taxID, row[5], genomeSize, potentialCoverage)
            results[kr.taxID] = kr
    return results

###########################
### Kraken Results data ###
###########################
class KrakenResult:
    def __init__(self, p, r1, r2, rc, tax, name, gs, cov):
        self.percentageOfReadsMappingToThis = float(p)
        self.readDepthThisLevelOnly = int(r1)
        self.readDepthThisOrLowerLevel = int(r2)
        # D superkingdom, K kingdom, P phylum, C class, O order, F family, G genus, S species see https://github.com/DerrickWood/kraken/blob/master/scripts/kraken-report
        self.rankCode = rc
        self.taxID = tax
        self.scientificName = name.lstrip()
        self.genomeSize = gs
        # assume each read is 100bp, and that none of the reads are overalapping, what coverage would we then have?
        self.potentialCoverage = cov

    # this function will determine if the given kraken result passes the filters
    # filters include minimum read depth, whether or not its a pathogen,
    # which type of rank code (only genus, species, or no rank),
    # or the maximum potential coverage (under the assumptions)
    def passFilters(self, minThisReadDepth, minTotalReadDepth, minGenusReadDepth, whiteList, minPotentialCoverage):
        passedXFilters = 0
        failReason = ""
        genome_size = self.genomeSize
        if(re.match("G", self.rankCode)):
            if(self.readDepthThisLevelOnly>=minGenusReadDepth and self.readDepthThisOrLowerLevel>=minTotalReadDepth):
                passedXFilters = passedXFilters + 1
            else:
                failReason = failReason + "NotEnoughReadDepth; "
        else:
            if(self.readDepthThisLevelOnly>=minThisReadDepth and self.readDepthThisOrLowerLevel>=minTotalReadDepth):
                passedXFilters = passedXFilters + 1
            else:
                failReason = failReason + "NotEnoughReadDepth; "
        if(self.taxID in whiteList.isWhitelisted):
            passedXFilters = passedXFilters + 1
        else:
            failReason = failReason + "NotWhiteListed; "
        # Make sure to not include species that don't have a genome size
        # Instead include 'no rank' taxa that have a genome size
        if self.rankCode in ["S","G"]:
            if self.rankCode == "S":
                if genome_size:
                    passedXFilters = passedXFilters + 1
                else:
                    failReason = failReason + "NoGenomeSizeForSpecies;"
            else:
                passedXFilters = passedXFilters + 1
        else:
            if genome_size:
                passedXFilters = passedXFilters + 1
                # If it has a genome size, then for all practical purposes it can be considered a species
                self.rankCode = "S"
            else:
                failReason = failReason + "NotASpeciesOrGenus;"
        if(self.potentialCoverage == "" or self.potentialCoverage>=minPotentialCoverage):
            passedXFilters = passedXFilters + 1
        else:
            failReason = failReason + "NotEnoughCoverage; "
        if(passedXFilters>=4):
            return (True, "")
        return (False, failReason);

    def printAsTabDelLine(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.taxID,self.readDepthThisLevelOnly,self.readDepthThisOrLowerLevel,self.scientificName,self.genomeSize,self.potentialCoverage,self.rankCode)

######################################
### The raw output on a read level ###
######################################
class CentrifugeReadLevelClassifications:
    def __init__(self, fName, printStat, whitelisted, tree):
        self.species = {}
        self.readStat = {}
        foundTaxIDs = {}
        readLengths = {}
        checked_taxa = {}
        with open(fName, 'r') as fileIn:
            for index, line in enumerate(fileIn, start=1):
                if index == 1:
                    continue
                if index%1000 == 0:
                    logging.info("{} lines read. {} reads parsed ({} taxID checked)".format(index,len(readLengths),len(checked_taxa)))
                row = line.strip().split("\t")
                readID, seqID, taxID, queryLength = row[0], row[1], row[2], row[6]
                # If taxID has been parsed already, get list of reads assigned
                # If not, start new list
                if taxID in self.species:
                    l = self.species[taxID]
                else:
                    l = []
                l.append(readID)

                # Perform a check for this taxID to see if it is whitelisted or has a child taxa whitelisted

                # Check if taxID has been retrieved already
                try:
                    haveAChildWhitelisted = checked_taxa[taxID]
                except KeyError:
                    # check if I have a child that is on the whitelist
                    try:
                        children = tree.get_descendant_taxa(int(taxID))
                    except ValueError:
                        checked_taxa[taxID] = False
                        continue

                    if len(set(whitelisted).intersection(children)) > 0:
                        haveAChildWhitelisted = True
                    else:
                        haveAChildWhitelisted = False
                    checked_taxa[taxID] = haveAChildWhitelisted

                if printStat or (taxID in whitelisted) or haveAChildWhitelisted:
                    self.species[taxID] = l

                # Store stats for read
                try:
                    foundTaxIDs[readID].append(taxID)
                except KeyError:
                    foundTaxIDs[readID] = [taxID]
                try:
                    readLengths[readID].append(int(queryLength))
                except KeyError:
                    readLengths[readID] = [int(queryLength)]

        logging.info("{} lines read. {} reads parsed ({} taxID checked)".format(index, len(readLengths),
                                                                                        len(checked_taxa)))
        # Collate stats for reads
        for readID in foundTaxIDs.keys():
            mean_len = int(numpy.round(numpy.mean(readLengths[readID])))
            taxa_found = len(set(foundTaxIDs[readID]))
            self.readStat[readID] = (mean_len, taxa_found)

    def getStatistics(self, taxID):
        if taxID in self.species:
            # Get all reads assigned to this species
            listOfReadsForThisTax = self.species[taxID]
            numberOfDifferentSpecies = []
            readLengths = []
            for readID in listOfReadsForThisTax:
                temp = self.readStat[readID]
                numberOfDifferentSpecies.append(temp[1])
                readLengths.append(temp[0])
            twentyFive = numpy.percentile(numberOfDifferentSpecies, 25)
            median = numpy.percentile(numberOfDifferentSpecies, 50) # MedianSpeciesPerRead
            seventyFive = numpy.percentile(numberOfDifferentSpecies, 75)
            med = numpy.percentile(readLengths, 50) # MedianReadLengthForSpecies
            stats = [twentyFive, median, seventyFive, med]
            return stats
        else:
            return [0,0,0,0]

class KrakenReadLevelClassifications:
    def __init__(self, fName, printStat, whitelisted, tree):
        # for each species, make a list of the reads assigned to it
        self.species = {}
        # for each read, remember both the number of species, and length
        self.readStat = {}
        checked_taxa = {}
        with open(fName, 'r') as fileIn:
            #reader = csv.reader(fileIn, delimiter='\t')
            for index, line in enumerate(fileIn):
                row = line.strip().split("\t") # remove the newline and split on .tab
            #for row in reader:
                if row[0] == "U":
                    continue # if the result is unclassified, dont bother interpreting it
                readID = row[1]

                lcaMappedKMerAssignments = row[5] # for example "562:13 561:4 A:31 0:1 562:3"
                setOfAssignments = lcaMappedKMerAssignments.split(" ") # for example 562:13
                #if(printStat):
                #    self.readLengths[readID] = int(row[3])
                foundTaxIDs = {}
                for a in setOfAssignments:
                    m = a.split(':')
                    taxID = m[0]
                    if (taxID == 0) | (str(taxID) == "0"):
                        continue # this k-mer was not in the database
                    if taxID == "A":
                        continue # this k-mer contained an ambigous nucleotide
                    l = [] # make an empty list
                    if taxID in self.species:
                        l = self.species[taxID]
                    l.append(readID)
                    try:
                        haveAChildWhitelisted = checked_taxa[taxID]
                    except KeyError:
                        # check if I have a child that is on the whitelist
                        try:
                            children = tree.get_descendant_taxa(int(taxID))
                        except ValueError:
                            checked_taxa[taxID] = False
                            continue
                        if len(set(whitelisted).intersection(children)) > 0:
                            haveAChildWhitelisted = True
                        else:
                            haveAChildWhitelisted = False
                        checked_taxa[taxID] = haveAChildWhitelisted
                    if printStat or (taxID in whitelisted) or haveAChildWhitelisted:
                        self.species[taxID] = l
                    # only remember the species if I should print statistics, or if the current one is a pathogen, or if its parent is a pathogen
                    if printStat or (taxID in whitelisted) or haveAChildWhitelisted:
                        self.species[taxID] = l
                    foundTaxIDs[taxID] = 1
                #if (printStat):
                self.readStat[readID] = (int(row[3]), len(foundTaxIDs))

    def getStatistics(self, taxID):
        if(taxID in self.species):
            # Get all reads assigned to this species
            listOfReadsForThisTax = self.species[taxID]
            numberOfDifferentSpecies = []
            readLengths = []
            for readID in listOfReadsForThisTax:
                temp = self.readStat[readID]
                numberOfDifferentSpecies.append(temp[1])
                readLengths.append(temp[0])
            a = numpy.array(numberOfDifferentSpecies)
            twentyFive = numpy.percentile(a, 25)
            median = numpy.percentile(a, 50) # MedianSpeciesPerRead
            seventyFive = numpy.percentile(a, 75)
            b = numpy.array(readLengths)
            med = numpy.percentile(b, 50) # MedianReadLengthForSpecies
            stats = [twentyFive, median, seventyFive, med]
            return(stats)
        else:
            return([0,0,0,0])


################################################################
### keep tracks of only NCBI taxonomy id and scientific name ###
################################################################
class WhiteList:
    def __init__(self, whiteListFile):
        self.isWhitelisted = {}
        with open(whiteListFile, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                taxID = row[0]
                if len(row) >1:
                    taxName = row[1]
                    self.isWhitelisted[taxID] = taxName
                else:
                    self.isWhitelisted[taxID] = taxID

#########################################################
### for each genome from NCBI there is an genome size ###
### of course there are none for the higher levels!   ###
#########################################################
class GenomeSizes:
    def __init__(self, genomeSizeFile):
        self.genomeSizes = {}
        with open(genomeSizeFile, 'r') as f:
            reader = csv.reader(f, delimiter = '\t')
            for row in reader:
                self.genomeSizes[row[0]] = int(row[1])

    def get(self, taxID):
        # there are genomes that have no genome size in these NCBI reports
        if taxID in self.genomeSizes:
            return self.genomeSizes[taxID]
        return ""

    def coverage(self, taxID, nrReads):
        # for the genera it's not really possible to calculate coverage
        if taxID in self.genomeSizes:
            return (nrReads*100) / self.genomeSizes[taxID]
        return ""