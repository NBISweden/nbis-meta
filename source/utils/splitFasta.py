import re, logging, os
from argparse import ArgumentParser
from Bio.SeqIO import parse
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s:%(levelname)s:%(message)s')


def store_seqs(f):
    seqs_taxid = {}
    seqs_acc = {}
    for record in parse(f, 'fasta'):
        try:
            id = (record.id).split("|")[1]
        except IndexError:
            # Likely the updated accession
            id = (record.id).split(".")[0]
        rec = ">{}\n{}\n".format(record.description, str(record.seq))
        if (re.match("kraken:taxid", record.id)):
            try:
                seqs_taxid[id].append(rec)
            except KeyError:
                seqs_taxid[id] = [rec]
        else:
            seqs_acc[id] = rec
    return seqs_acc, seqs_taxid


def id2taxid(ids, mapfiles=[]):
    idtaxid = {}
    total = len(ids)
    found = 0
    if len(mapfiles)==0:
        logging.info("No mapfile, attempting to map ids to taxids using Entrez efetch...")
        import urllib3
        urllib3.disable_warnings()
        http = urllib3.PoolManager()

        for id in ids:
            s = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=xml".format(id)
            r = http.request('GET', s)
            m = re.search("TSeq_taxid>\d+<\/TSeq_taxid>", str(r.data))
            try:
                taxid = m.group().split("<")[0].split(">")[-1]
            except AttributeError:
                print(s)

            if taxid:
                found += 1
            else:
                taxid = "_missing"
            idtaxid[id] = taxid
            if found%100 == 0:
                logging.info("Found {}/{} taxids ({}% done)".format(found,total,round(float(found)/total*100),2))

    else:
        for mapfile in mapfiles:
            logging.info("Using {} as mapfile to map ids to taxids".format(mapfile))
            with open(mapfile, 'r') as fh:
                for line in fh:
                    items = (line.rstrip()).rsplit()
                    if len(items) == 2:
                        acc,taxid = items
                        acc = acc.split(".")[0]
                        id = acc_ver = ""
                    elif len(items) == 4:
                        acc,acc_ver,taxid,id = items
                    if acc in ids:
                        idtaxid[acc] = taxid
                    elif id in ids:
                        idtaxid[id] = taxid
                    else:
                        continue
                    idtaxid[id] = taxid
                    found += 1
                    if found%100 == 0:
                        logging.info("Found {}/{} taxids ({}% done)".format(found, total, round(float(found) / total * 100), 2))
                    if found == total:
                        break
    missing = set(ids).difference(set(idtaxid.keys()))
    logging.info("Found {}/{} taxids ({} missing)".format(found,total,len(missing)))
    for id in missing:
        idtaxid[id] = "_missing"
    return idtaxid


def add_deprecated(mapfile, idtaxid):
    logging.info("Updating from deprecated ids")
    updated = 0
    with open(mapfile, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            dep_id,id = line.rsplit()
            if id in idtaxid:
                idtaxid[dep_id] = idtaxid[id]
                updated+=1
    logging.info("{} ids updated".format(updated))
    return idtaxid

def write(seqs_taxid, outdir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    logging.info("Writing sequence records to files in {}/".format(outdir))
    for taxid, records in seqs_taxid.items():
        with open("{}/{}.fa".format(outdir,taxid), 'w') as fh:
            for record in records:
                fh.write(record)

def main():
    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True,
                        help="Fasta file to split")
    parser.add_argument("-m", "--mapfiles", nargs="+",
                        help="Files mapping accession numbers to tax ids")
    parser.add_argument("-d", "--deprecated",
                        help="File with deprecated gi mappings")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Directory to store split fasta files. Defaults to current directory.")
    args = parser.parse_args()
    seqs_acc, seqs_taxid = store_seqs(args.fasta)
    logging.info("Stored {} sequences with taxid, {} sequences with gi/accession".format(len(seqs_taxid),len(seqs_acc)))
    if len(seqs_acc) > 0:
        idtaxid = id2taxid(seqs_acc.keys(), args.mapfiles)
        if args.deprecated:
            idtaxid = add_deprecated(args.deprecated, idtaxid)
        for id, records in seqs_acc.items():
            try:
                seqs_taxid[idtaxid[id]] += records
            except KeyError:
                seqs_taxid[idtaxid[id]] = records
    write(seqs_taxid, args.outdir)


if __name__ == '__main__':
    main()