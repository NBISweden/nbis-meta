from ete3 import NCBITaxa
import os, pandas as pd, subprocess, sys, gzip as gz
import io
from argparse import ArgumentParser
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s:%(levelname)s:%(name)s:%(message)s'
)

def read_lines(fh):
    all_prot_mapped_taxids = {}
    for i, line in enumerate(fh):
        if i == 0:
            continue
        try:
            _, _, taxid, _ = line.decode().split("\t")
        except AttributeError:
            _, _, taxid, _ = line.split("\t")
        try:
            all_prot_mapped_taxids[taxid] += 1
        except KeyError:
            all_prot_mapped_taxids[taxid] = 1
    return all_prot_mapped_taxids


def get_all_taxaid(acc_to_prot_file):
    if ".gz" in acc_to_prot_file:
        with gz.open(acc_to_prot_file) as prot_map_fh:
            with io.BufferedReader(prot_map_fh) as prot_map_buff_fh:
                return read_lines(prot_map_buff_fh)
    else:
        with open(acc_to_prot_file) as prot_map_fh:
            return read_lines(prot_map_fh)


def get_superkingdom(lineage, ncbi_taxa):
    for t in lineage:
        t = int(t)
        rank = ncbi_taxa.get_rank([t])[t]
        name = ncbi_taxa.get_taxid_translator([t])[t]
        if rank == "superkingdom":
            return name
    return "Unknown"


def count_ranks(taxid_counts, ncbi_taxa):
    superkingdoms = ["Bacteria", "Eukaryota", "Archaea", "Viruses", "Unknown"]
    result = {}
    for s in superkingdoms:
        result[s] = {"phylum": set(), "class": set(), "order": set(), "family": set(), "genus": set(), "species": set()}

    for taxid in taxid_counts.keys():
        try:
            lineage = ncbi_taxa.get_lineage(taxid)
        except ValueError as err:
            logging.info("Warning: {e}".format(e=err))
            continue

        superkingdom = get_superkingdom(lineage, ncbi_taxa)
        for t in lineage:
            t = int(t)
            rank = ncbi_taxa.get_rank([t])[t]
            name = ncbi_taxa.get_taxid_translator([t])[t]
            try:
                result[superkingdom][rank].add(name)
            except KeyError:
                continue
    for s in superkingdoms:
        for rank in result[s].keys():
            result[s][rank] = len(result[s][rank])
    return result

def main():
    parser = ArgumentParser()
    parser.add_argument("-m", "--map", help="Protein accession -> taxid map file", nargs="+")
    parser.add_argument("-d", "--dbfile",
        help="Path to sqlite database to create using ete3. If not specified it will be \
        saved as {home}/.etetoolkit/taxa.sqlite".format(home=os.path.expanduser('~')))

    args = parser.parse_args()

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

    results = {}
    for i,acc_to_prot_file in enumerate(args.map):
        logging.info("Getting taxaids from {f}".format(f=acc_to_prot_file))
        fname = os.path.basename(acc_to_prot_file)
        taxid_counts = get_all_taxaid(acc_to_prot_file)
        logging.info("Counting ranks...")
        rank_counts = count_ranks(taxid_counts, ncbi_taxa)
        results[fname] = rank_counts.copy()
        logging.info("done.")

    df = pd.DataFrame()
    for fname in results.keys():
        for superkingdom in sorted(results[fname].keys()):
            df_ = pd.DataFrame(results[fname][superkingdom], index=[superkingdom])
            df_.index.name = "superkingdom"
            df_ = df_.reset_index()
            df_.index = [fname]
            df = pd.concat([df, df_])
    df = df[["superkingdom","phylum","class","order","family","genus","species"]]
    df.index.name = "database"
    df.to_csv(sys.stdout, sep="\t")

if __name__ == '__main__':
    main()

