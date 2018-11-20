# coding: utf-8
"""
A script to create a tab seperated file containing the taxonomy for all leaves in the ncbi taxonomy
defined for the rank levels superkingdom, phylum, class, order, family, genus and species.

Modified from https://github.com/EnvGen/toolbox/blob/master/scripts/assign_taxonomy_from_blast/new_lineage_file.py
"""
from ete3 import NCBITaxa
import argparse
import sys
import gzip as gz
import logging
import subprocess
import os
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s:%(levelname)s:%(name)s:%(message)s'
)

import io

def read_lines(fh):
    all_prot_mapped_taxids = set()
    for i,line in enumerate(fh):
        if i==0:
            continue
        try: _, _, taxid, _ = line.decode().split("\t")
        except AttributeError:  _, _, taxid, _ = line.split("\t")
        all_prot_mapped_taxids.add(int(taxid))
    return all_prot_mapped_taxids

def get_all_taxaid(acc_to_prot_file):
    logging.info("Fetching taxaids from accession to prot file")
    if ".gz" in acc_to_prot_file:
        with gz.open(acc_to_prot_file) as prot_map_fh:
            with io.BufferedReader(prot_map_fh) as prot_map_buff_fh:
                return read_lines(prot_map_buff_fh)
    else:
        with open(acc_to_prot_file) as prot_map_fh:
            return read_lines(prot_map_fh)

def find_base_annotation(last_annotated_level, i, fixed_lineage, names, full_lineage):
    base_annotation = None
    if last_annotated_level != -1:
        for j in  range(i,last_annotated_level, -1):
            if fixed_lineage[last_annotated_level] == names[full_lineage[j]]:
                if base_annotation is None:
                    base_annotation = fixed_lineage[last_annotated_level]
                return base_annotation
            else:
                base_annotation = names[full_lineage[j]]
    else:
        base_annotation = names[full_lineage[0]]
    if base_annotation == "root":
        # If the full lineage consists of only root, use Unclassified as base_annotation
        # This can occur if taxids are taken from a clustered database such as Uniref90
        if len(full_lineage) == 1:
            base_annotation = 'Unclassified'
        else:
            base_annotation = names[full_lineage[1]]

    return base_annotation

def add_annotation(last_annotated_level, i, fixed_lineage, names, full_lineage, lin_taxaid, missed_levels):
    base_annotation = find_base_annotation(last_annotated_level, i, fixed_lineage, names, full_lineage)

    for level in range(last_annotated_level+1, last_annotated_level+missed_levels):
        fixed_lineage.append("{}".format(base_annotation))
    fixed_lineage.append(names[lin_taxaid])
    last_annotated_level += missed_levels
    return fixed_lineage, last_annotated_level

def taxaid_to_fixed_lineage(taxaid, ncbi):
    fixed_ranks = {"superkingdom": 0,
                  "phylum": 1,
                  "class": 2,
                  "order": 3,
                  "family": 4,
                  "genus": 5,
                  "species": 6}

    try:
        full_lineage = ncbi.get_lineage(taxaid)
    except ValueError:
        return []
    if not full_lineage:
        return []
    names = ncbi.get_taxid_translator(full_lineage)
    ranks = ncbi.get_rank(full_lineage)

    fixed_lineage = []
    last_annotated_level = -1

    for i, lin_taxaid in enumerate(full_lineage):

        if ranks[lin_taxaid] in fixed_ranks:
            if fixed_ranks[ranks[lin_taxaid]] == last_annotated_level:
                # Dealing with cases such as taxaid: 1918211
                # which have two values for the same rank (genus)
                continue
            missed_levels = fixed_ranks[ranks[lin_taxaid]] - last_annotated_level

            if missed_levels == 1: # No missed level
                fixed_lineage.append(names[lin_taxaid])
                last_annotated_level += 1
            elif last_annotated_level == -1: # Missing value for superkingdom:
                fixed_lineage, last_annotated_level = add_annotation(last_annotated_level, i, fixed_lineage, names, full_lineage, lin_taxaid, missed_levels)
            else:
                # Are there any items in the lineage that are not on a fixed level
                fixed_lineage, last_annotated_level = add_annotation(last_annotated_level, i, fixed_lineage, names, full_lineage, lin_taxaid, missed_levels)

    if last_annotated_level != 6: # Species level is missing in lineage
        base_annotation = find_base_annotation(last_annotated_level, i, fixed_lineage, names, full_lineage)

        for level in range(last_annotated_level+1, 7):
            fixed_lineage.append("{}".format(base_annotation))

    if len(fixed_lineage) != 7:
        logging.warning("Lineage for {} is of length {}".format(taxaid, len(fixed_lineage)))

    assert len(fixed_lineage) == 7

    return fixed_lineage

def write_lineage(taxaids, ncbi_taxa, outfile):
    if not outfile:
        fhout = sys.stdout
    else:
        fhout = open(outfile, 'w')
    logging.info("Print lineage for each taxaid")
    for taxaid in taxaids:
        lineage = taxaid_to_fixed_lineage(taxaid, ncbi_taxa)
        if lineage:
            fhout.write(str(taxaid) + "\t")
            fhout.write("\t".join(lineage)+ "\n")
        else:
            logging.warning("Taxaid missing in ncbi taxonomy: {}".format(taxaid))

def main(args):
    logging.info("Loading ncbi taxonomy")
    if not args.dbfile:
        ncbi_taxa = NCBITaxa()
    else:
        dirname = os.path.dirname(args.dbfile)
        cmd = "mkdir -p "+dirname
        p1 = subprocess.Popen(cmd, shell=True, stdin=None)
        p1.wait()
        cmd = "touch "+args.dbfile
        p2 = subprocess.Popen(cmd, shell=True, stdin=None)
        p2.wait()
        ncbi_taxa = NCBITaxa(args.dbfile)

    taxaids = get_all_taxaid(args.accession_to_taxid)

    write_lineage(taxaids, ncbi_taxa, args.outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-a", "--accession_to_taxid",
                        help="The accession to protein file prot.accession2taxid file from ncbi: http://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/")
    parser.add_argument("-o", "--outfile",
                        help="Write lineage to outfile. Defaults to stdout.")
    parser.add_argument("-d", "--dbfile",
                        help="Path to sqlite database to create using ete3. If not specified it will be saved as ~/.etetoolkit/taxa.sqlite")

    args = parser.parse_args()

    main(args)
