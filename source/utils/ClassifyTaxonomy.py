import argparse
import logging
import operator
from collections import Counter
from collections import defaultdict

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def calculate_prot_length(start, end):
    return (abs(int(end) - int(start)) + 1)/3.0


def test_calculate_prot_length():
    assert calculate_prot_length(12, 1) == 4
    assert calculate_prot_length(1, 12) == 4


def calculate_fHit(percIdentity, queryEnd, queryStart, qLength):
    '''Calculates percent identity normalized to the fraction of the alignment to the query'''
    alnLength_in_query = abs(int(queryEnd) - int(queryStart)) + 1
    fHit = float(alnLength_in_query) / qLength
    fHit *= float(percIdentity) / 100.0
    fHit = min(1.0, fHit)
    return fHit


def test_fHit_calculation():
    percIdentity = 90.0
    queryEnd = 309
    queryStart = 10
    qLength = 600
    ## 309-10+1 = 300, 300/600 = 0.5, 0.5*(90/100) = 0.9*0.5 = 0.45
    assert calculate_fHit(percIdentity, queryEnd, queryStart, qLength) == 0.45


def read_blast_lines(fh, lengths, min_id):
    # k191_83_2 gi|973180054|gb|KUL19018.1| 71.2 73 21 0 9 81 337 409 6.6e-24 118.2
    # qureryId subjectId percIdentity alnLength mismatchCount gapOpenCount queryStart queryEnd subjectStart subjectEnd eVal bitScore

    matches = defaultdict(list)
    gids = Counter()
    for line in fh:
        line = line.rstrip()
        (queryId, subjectId, percIdentity, alnLength, mismatchCount,
         gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore) = line.split("\t")

        gid = subjectId
        ## Get length of query
        qLength = lengths[queryId]
        ## Calculate percent identity normalized to alignment length
        fHit = calculate_fHit(percIdentity, queryEnd, queryStart, qLength)

        if float(percIdentity) >= min_id:
            matches[queryId].append((gid, fHit))
            gids[gid] += 1
    return (matches, gids)


def read_blast_input(blastinputfile, lengths, min_id):
    with open(blastinputfile) as fh:
        (matches, gids) = read_blast_lines(fh, lengths, min_id)
    return (matches, gids.keys())


def run_map_back(mapBack, tokens):
    for depth in range(6, 0, -1):
        if tokens[depth] not in mapBack[depth] and tokens[depth] != 'None':
            for depth2 in range(depth - 1, -1, -1):
                mapBack[depth][tokens[depth]].append(tokens[depth2])
    return mapBack


def read_lineage_lines(fh):
    mapping = {}
    mapBack = defaultdict(lambda: defaultdict(list))
    for line in fh:
        line = line.rstrip()
        tokens = line.split("\t")
        (taxaid, domainid, phylumid, classid, orderid, familyid, genusid, speciesid) = tokens
        ## Map taxaid to the lineage
        mapping[int(taxaid)] = [domainid, phylumid, classid, orderid, familyid, genusid, speciesid];
        tokens.pop(0)  ## Remove taxaid from tokens
        ## Iterate over the ranks and map each rank to the full previous lineage
        mapBack = run_map_back(mapBack, tokens)
    return (mapping, mapBack)


def read_lineage_file(lineage_file):
    with open(lineage_file) as fh:
        (mapping, mapBack) = read_lineage_lines(fh)
    return (mapping, mapBack)


def read_query_length_lines(fh):
    lengths = {}
    for line in fh:
        line = line.rstrip()
        (queryid, length) = line.split("\t")
        lengths[queryid] = float(length)
    return lengths


def read_query_length_file(query_length_file):
    with open(query_length_file) as fh:
        lengths = read_query_length_lines(fh)
    return lengths


def map_accessions(accs, fh):
    total_accs = len(accs)
    found = 0
    mappings = dict([(acc, -1) for acc in accs])
    for i, line in enumerate(fh):
        if i == 0:
            continue
        _, acc_ver, taxid, _ = line.split("\t")
        # Only add taxids for the given acc
        if acc_ver in mappings:
            mappings[acc_ver] = int(taxid)
            found+=1
            if found==total_accs:
                break
    return mappings


def read_accessions_file(accs, mapping_file):
    with open(mapping_file) as fh:
        mappings = map_accessions(accs, fh)
    return mappings


def read_bed_lines(fh):
    contigGenes = defaultdict(list)
    contigLengths = Counter()
    lengths = {}
    for line in fh:
        line = line.rstrip()
        items = line.split("\t")
        # GFF file in the format used for htseq
        if len(items) == 9:
            contig,start,end = items[0],items[3],items[4]
            gene_id = items[-1]
            if " " in gene_id:
                gene = contig+"_"+gene_id.split("_")[-1]
        elif len(items) == 4:
            contig, start, end, gene = line.split("\t")
        contigGenes[contig].append(gene)
        length = calculate_prot_length(start, end)
        lengths[gene] = length
        contigLengths[contig] += length
    return (lengths, contigGenes, contigLengths)


def read_bed_file(bedfile):
    with open(bedfile) as fh:
        (lengths, contigGenes, contigLengths) = read_bed_lines(fh)
    return lengths, contigGenes, contigLengths


def calculate_taxa_weight(fHit, min_id_taxa):
    ## If min_id_taxa = 0.95, fHit needs to be > 0.95 for this depth
    weight = (fHit - min_id_taxa) / (1.0 - min_id_taxa)
    weight = max(weight, 0.0)
    return weight


def test_calculate_taxa_weight():
    fHit = 0.75
    min_id_taxa = 0.5
    assert calculate_taxa_weight(fHit, min_id_taxa) == 0.5


def collate_gene_hits(matchs, mapping, lineages, taxa_identities):
    collate_hits = list()
    for depth in range(7):
        collate_hits.append(Counter())

    added_matches = set()
    missing_taxaid = {}

    ## Iterate the hits (protein accessions) and normalized identity, sorted by decreasing normalized percent identity
    for (match, fHit) in sorted(matchs, key=lambda x: x[1], reverse=True):
        if mapping[match] <= -1:
            continue
        ## Get the taxid and lineage for the protein accession
        tax_id = mapping[match]

        ## Only add best hit per species
        if tax_id in added_matches:
            continue
        added_matches.add(tax_id)

        if tax_id not in lineages and tax_id not in missing_taxaid.keys():
            logging.warning("Taxa id {} is missing from lineage file".format(tax_id))
            missing_taxaid[tax_id] = ""
            continue

        hits = lineages[tax_id]
        for depth in range(7):
            if hits[depth] != "None":
                ## Calculate the normalized weight for each depth, if < 0 the hit is not included
                weight = calculate_taxa_weight(fHit, taxa_identities[depth])
                if weight > 0.0:
                    collate_hits[depth][hits[depth]] += weight  # could put a transform in here
    return collate_hits


def make_assignment(geneAssign, gene, collate_hits, mapBack, min_fraction):
    for depth in range(6, -1, -1):
        collate = collate_hits[depth]
        dWeight = sum(collate.values())
        sortCollate = sorted(collate.items(), key=operator.itemgetter(1), reverse=True)
        nL = len(collate)  ## Number of hits collated at depth
        if nL > 0:
            dP = 0.0
            ## dWeight is the sum of weights at this depth
            if dWeight > 0.0:
                ## dP is the weight of the best hit at depth, divided by the total weight at depth
                dP = float(sortCollate[0][1]) / dWeight
                ## If the dP for the best hit is larger than min_fraction then assign at this depth and map back
                if dP > min_fraction:
                    geneAssign[gene][depth] = (sortCollate[0][0], dP)
                    ## Get back lineage for the best assignment
                    assignBack = mapBack[depth][sortCollate[0][0]]
                    depth2 = depth - 1
                    for assignB in assignBack:
                        geneAssign[gene][depth2] = (assignB, 1.0)
                        depth2 -= 1
                    break

                else:
                    geneAssign[gene][depth] = ('Unclassified', -1.0)
            ## If dWeight is 0.0, set 'Unclassified'
            else:
                geneAssign[gene][depth] = ('Unclassified', -1.0)
        ## If no hits collated at depth, set 'Unclassified'
        else:
            geneAssign[gene][depth] = ('Unclassified', -1.0)
    return geneAssign


def assign_unclassified():
    d = {}
    for depth in range(7):
        d[depth] = ('Unclassified',-1.0)
    return d


def assign_taxonomy(matches, mapping, mapBack, lineages, lengths, contigGenes, min_fraction, taxa_identities):
    geneAssign = defaultdict(dict)
    geneassign_from_contig = defaultdict(dict)
    contigAssign = defaultdict(dict)
    total_genes = len(lengths)
    assigned_genes = 0
    for contig, genes in contigGenes.items():
        collate_hits = list()
        for depth in range(7):
            collate_hits.append(Counter())
        for gene in genes:
            assigned_genes+=1
            if assigned_genes % 100000 == 0:
                assigned_percent = round(float(assigned_genes)/total_genes*100)
                logging.info("{}/{} assigned genes ({}%)".format(assigned_genes,total_genes,assigned_percent))
            ## Perform taxonomic assignments for all genes on contig
            try:
                matchs = matches[gene]
                collated_gene_hits = collate_gene_hits(matchs, mapping, lineages, taxa_identities)
                geneAssign = make_assignment(geneAssign, gene, collated_gene_hits, mapBack, min_fraction)
            except KeyError:
                logging.info("No matches found for "+gene)
                geneAssign[gene] = assign_unclassified()
                continue
            ## Then add assignments to the contig for collating
            for depth in range(7):
                try:
                    (assignhit, genef) = geneAssign[gene][depth]
                except KeyError:
                    continue
                if assignhit != 'Unclassified':
                    collate_hits[depth][assignhit] += lengths[gene]

        # Contigs are assigned a taxonomy based on LCA for
        # all genes assigned on at least kingdom level.
        dWeight = sum(collate_hits[0].values())
        for depth in range(7):
            collate = collate_hits[depth]
            sortCollate = sorted(collate.items(), key=operator.itemgetter(1), reverse=True)
            nL = len(collate)
            if nL > 0:
                dP = 0.0
                if dWeight > 0.0:
                    dP = float(sortCollate[0][1]) / dWeight
                    if dP > min_fraction:
                        contigAssign[contig][depth] = (sortCollate[0][0], dP)#, sortCollate[0][1])
                        for gene in genes: geneassign_from_contig[gene][depth] = (sortCollate[0][0], dP)
                    else:
                        contigAssign[contig][depth] = ('Unclassified', 0.)
                        for gene in genes: geneassign_from_contig[gene][depth] = ('Unclassified', 0.)
                else:
                    contigAssign[contig][depth] = ('Unclassified', 0.)
                    for gene in genes: geneassign_from_contig[gene][depth] = ('Unclassified', 0.)
            else:
                contigAssign[contig][depth] = ('Unclassified', 0.)
                for gene in genes: geneassign_from_contig[gene][depth] = ('Unclassified', 0.)
    return contigAssign, geneAssign, geneassign_from_contig


def write_assigns(assign, assign_file, support_file):
    with open(assign_file, "w") as assign_fh, open(support_file, "w") as support_fh:
        for key in assign.keys():
            assign_fh.write('%s' % key)
            support_fh.write('%s' % key)
            last_known = ""
            for depth in range(7):
                (a, p) = assign[key][depth]
                if a!="Unclassified":
                    last_known = a
                else:
                    if last_known:
                        # Don't add "Unclassified" prefix to unclassified taxa names
                        if last_known.lower()[0:12] != "unclassified":
                            a = a + "." + last_known
                        else:
                            a = last_known
                assign_fh.write('\t%s' % a)
                support_fh.write('\t%.3f' % p)
            assign_fh.write('\n')
            support_fh.write('\n')
            assign_fh.flush()
            support_fh.flush()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--blast_input_file",
                        help="Blast input file with blast 6 matches to taxaid database *.b6")
    parser.add_argument('-a', '--acc_taxaid_mapping_file',
                        help="mapping from accession to taxaid gzipped")
    parser.add_argument('-b', '--bed_file', required=True,
                        help="bed file defining genes on contigs")
    parser.add_argument('-l', '--lineage_file', required=True,
                        help="text taxaid to lineage mapping")
    parser.add_argument('-f', '--min_fraction', default=0.5, type=float,
                        help="Minimum fraction of weights needed to assign to a particular level, default 0.5. 0.9 \
                        would be more strict.")
    parser.add_argument('-m', '--min_id', default=40.0, type=float,
                        help="Minimum allowed percent identity to parse a hit")
    parser.add_argument('-o', '--output_name', type=str, default="output",
                        help="String specifying output base name. Gene and contig assignments will be written \
                        with suffixes '_genes.tab and '_contigs.tab'. Support values are written with the additional \
                        suffixes _genes.supports.tab and _contigs.supports.tab")
    parser.add_argument("-t", "--taxa_identities", nargs="*", default=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95],
                        help="Specify normalized weights for the ranks: superkingdom phylum class order\
                         family genus species. Defaults to 0.4 0.5 0.6 0.7 0.8 0.9 0.95")

    args = parser.parse_args()

    taxa_identities = [float(x) for x in args.taxa_identities]

    if len(args.taxa_identities) != 7:
        raise Exception(
            "Please specify taxa identities as an array of length 7. Length of {} is {}.".format(args.taxa_identities, len(args.taxa_identities))
            )

    logging.info("Reading bed file")
    (lengths, contigGenes, contigLengths) = read_bed_file(args.bed_file)
    logging.info("Finished reading bed file")

    logging.info("Reading blast results file")
    (matches, gids) = read_blast_input(args.blast_input_file, lengths, args.min_id)
    logging.info("Finished reading in blast results file")

    logging.info("Reading lineage file")
    (lineages, mapBack) = read_lineage_file(args.lineage_file)
    logging.info("Finished reading in lineage file")

    logging.info("Reading taxaid map file")
    mapping = read_accessions_file(gids, args.acc_taxaid_mapping_file)
    logging.info("Finished loading taxaid map file")

    logging.info("Assigning taxonomy...")
    contigassign, geneassign, geneassign_from_contigs = assign_taxonomy(matches, mapping, mapBack, lineages, lengths, contigGenes, args.min_fraction, taxa_identities)
    logging.info("Finished assigning taxonomy")

    write_assigns(geneassign, args.output_name+"_genes.tab", args.output_name+"_genes.supports.tab")
    write_assigns(contigassign, args.output_name + "_contigs.tab", args.output_name + "_contigs.supports.tab")
    write_assigns(geneassign_from_contigs, args.output_name + "_genes_from_contigs.tab", args.output_name\
                  + "_genes_from_contigs.supports.tab")
    logging.info("Results written to "+args.output_name+"_genes.tab "+ args.output_name+"_contigs.tab"\
                 +args.output_name+"_genes_from_contigs.tab")
    logging.info("Supports written to "+args.output_name+"_genes.supports.tab "+ args.output_name+\
                 "_contigs.supports.tab"+ args.output_name+"_genes_from_contigs.supports.tab")

if __name__ == "__main__":
    main()
