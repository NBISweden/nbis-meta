#!/usr/bin/env python
import logging
from time import time
import gzip as gz
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def calculate_norm_id(percIdentity, queryEnd, queryStart, qLength):
    '''Calculates percent identity normalized to the fraction of the alignment to the query'''
    alnLength_in_query = abs(int(queryEnd) - int(queryStart)) + 1
    norm_id = float(alnLength_in_query) / int(qLength)
    norm_id *= float(percIdentity) / 100.0
    norm_id = min(1.0, norm_id)
    return norm_id


def get_ranks(taxid, ncbitax):
    ranks = ["superkingdom","phylum","class","order","family","genus","species"]
    rank_dict = {}
    try:
        lineage = ncbitax.get_lineage(taxid)
    except ValueError:
        return False
    r = ncbitax.get_rank(lineage)
    for id, rank in r.items():
        if rank in ranks:
            rank_dict[id] = rank
    return rank_dict


def parse_query(query_res, query_taxid, ncbitax, norm_id_dist, rank_dict):
    try:
        query_ranks = rank_dict[query_taxid]
    except KeyError:
        query_ranks = get_ranks(query_taxid, ncbitax)
        rank_dict[query_taxid] = query_ranks
    for line in query_res:
        qseqid, sseqid, pident, qlen, qstart, qend, length, bitscore, evalue, staxids = line
        # Get the normalized percent id
        norm_id = calculate_norm_id(pident, qend, qstart, qlen)
        try:
            subject_ranks = rank_dict[staxids]
        except KeyError:
            subject_ranks = get_ranks(staxids, ncbitax)
            rank_dict[staxids] = subject_ranks
        if not subject_ranks:
            continue
        # Get the difference between the query ranks and the subject ranks
        diff = set(query_ranks.keys()).difference(set(subject_ranks.keys()))
        rank_diffs = [query_ranks[taxid] for taxid in list(diff)]
        if len(rank_diffs) > 0:
            for rank in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
                # Store the normalized id at the highest rank that differs
                if rank in rank_diffs:
                    norm_id_dist[rank].append(norm_id)
                    break
    return norm_id_dist, rank_dict


def read_blast_file(f):
    from ete3 import NCBITaxa
    ncbitax = NCBITaxa()
    norm_id_dist = {"superkingdom": [], "phylum": [], "class": [], "order": [], "family": [], "genus": [], "species": []}
    query_res = []
    rank_dict = {}
    queries_parsed = 0
    current_query = ""
    query_taxid = ""
    with gz.open(f, 'rt') as fh:
        for line_num, line in enumerate(fh, start=1):
            # Get the values from the line
            qseqid, sseqid, pident, qlen, qstart, qend, length, bitscore, evalue, staxids = line.rstrip().rsplit()
            # Strip the UniRef100_ prefix
            qseqid = qseqid.replace("UniRef100_","")
            # If a new query is detected, parse the stored results
            if qseqid != current_query:
                queries_parsed+=1
                if queries_parsed % 100 == 0 and current_query:
                    logging.info("{} queries parsed".format(queries_parsed))
                if not current_query:
                    current_query = qseqid
                else:
                    if len(query_res) > 0:
                        logging.info("Parsing {} hits for query {} ({} lines parsed)".format(len(query_res),current_query, line_num))
                        start = time()
                        norm_id_dist,rank_dict = parse_query(query_res, query_taxid, ncbitax, norm_id_dist, rank_dict)
                        end = time()
                        time_per_res = (end-start)/float(len(query_res))
                        logging.info("{} hits parsed in {} seconds ({} s/result)".format(len(query_res), (end-start), time_per_res))


                    current_query = qseqid
                    query_res = []
            # If a self-hit, get the taxid for the query but don't store results
            if qseqid == sseqid:
                query_taxid = staxids
                continue
            # Add the results
            query_res.append([qseqid, sseqid, pident, qlen, qstart, qend, length, bitscore, evalue, staxids])
        # Handle last stored results
        norm_id_dist, rank_dict = parse_query(query_res, query_taxid, ncbitax, norm_id_dist, rank_dict)

    return norm_id_dist


def write_dist(norm_id_dist):
    import sys
    with sys.stdout as fh:
        for rank, dist in norm_id_dist.items():
            for item in dist:
                fh.write("{}\t{}\n".format(rank,item))


def get_taxidmap(f):
    map = {}
    with gz.open(f, 'rt') as fh:
        for line in fh:
            qseqid, _, taxid, _ = line.rstrip().split("\t")
            map[qseqid] = int(taxid)
    return map


def main():

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile",
                        help="Diamond blast output")
    args = parser.parse_args()
    logging.info("Parsing blast file")
    norm_id_dist = read_blast_file(args.infile)
    logging.info("Writing results")
    write_dist(norm_id_dist)

if __name__ == '__main__':
    main()