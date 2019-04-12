#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import json
import os
import urllib
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


# Download functions
####################

def get_kegg_module_hierarchy(s):
    hier = {}
    # First level is 'ko00002'
    for d1 in s['children']:
        c1 = d1['name']
        for d2 in d1['children']:
            c2 = d2['name']
            for d3 in d2['children']:
                c3 = d3['name']
                for module in d3['children']:
                    module_name = module['name']
                    kos = module['children']
                    hier[module_name] = {"Module_category1": c1, "Module_category2": c2,
                                         "Module_category3": c3, "KOs": []}
                    for ko in kos:
                        ko_name = ko['name']
                        ko = ko_name.split(" ")[0]
                        hier[module_name]["KOs"].append(ko)
    return hier


def get_kegg_ortholog_hierarchy(s):
    hier = {}
    # First level is 'ko00001'
    for d1 in s['children']:
        c1 = d1['name']
        for d2 in d1['children']:
            c2 = d2['name']
            for d3 in d2['children']:
                c3 = d3['name']
                if not "children" in d3.keys():
                    continue
                for ko in d3['children']:
                    ko_name = ko['name'].split("\t")[0]
                    ko_id = ko_name.split(" ")[0]
                    if "[EC:" in ko_name:
                        enzymes = ko_name.split("[")[-1].split("]")[0].lstrip("EC:").split(" ")
                    else:
                        enzymes = []
                    d = {"KO_category1": c1, "KO_category2": c2, "pathway": c3, "name": ko_name, "enzymes": enzymes}
                    try:
                        hier[ko_id].append(d)
                    except KeyError:
                        hier[ko_id] = [d]
    return hier


def get_kegg_ortholog_info(outdir):
    outdir = outdir.rstrip("/")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    url = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json"
    logging.info("Fetching ko00001.keg from www.kegg.jp")
    # Download file
    tmp_out = os.path.expandvars("$TMPDIR/ko00001.json")
    urllib.request.urlretrieve(url, tmp_out)
    with open(tmp_out) as fh:
        s = json.load(fh)
    hier = get_kegg_ortholog_hierarchy(s)
    pathways = {}
    ec2path = {}
    ko_out = "{}/kegg_kos.tsv".format(outdir)
    ko2path_out = "{}/kegg_ko2pathways.tsv".format(outdir)
    ko2ec_out = "{}/kegg_ko2ec.tsv".format(outdir)
    ec2path_out = "{}/kegg_ec2pathways.tsv".format(outdir)
    pathways_out = "{}/kegg_pathways.tsv".format(outdir)
    for f in [ko_out, ko2path_out, ko2ec_out, ec2path_out, pathways_out]:
        logging.info("Writing to {f}".format(f=f))
    # Write KEGG Ortholog names, KEGG Ortholog -> Pathway map, and KEGG Ortholog -> Enzyme map
    with open(ko_out, 'w') as fh_kos, open(ko2path_out, 'w') as fh_ko2path, open(ko2ec_out, 'w') as fh_ko2ec:
        fh_kos.write("ko\tKO_name\n")
        for ko_id, l in hier.items():
            for i, d in enumerate(l):
                if i == 0:
                    fh_kos.write("{}\t{}\n".format(ko_id, d["name"]))
                    for enzyme in d["enzymes"]:
                        fh_ko2ec.write("{}\t{}\n".format(ko_id, enzyme))
                fh_ko2path.write("{}\t{}\n".format(ko_id, "map{}".format(d["pathway"].split(" ")[0])))
                pathways[d["pathway"]] = {"Pathway_category1": d["KO_category1"],
                                          "Pathway_category2": d["KO_category2"]}
                for enzyme in d["enzymes"]:
                    try:
                        ec2path[enzyme].append("map{}".format(d["pathway"].split(" ")[0]))
                    except KeyError:
                        ec2path[enzyme] = ["map{}".format(d["pathway"].split(" ")[0])]
    # Write Pathway information
    with open(pathways_out, 'w') as fh_pathways:
        fh_pathways.write("{}\t{}\t{}\t{}\n".format("Pathway_id", "Pathway_name", "Pathway_category1",
                                                    "Pathway_category2"))
        for pathway, d in pathways.items():
            pathway_id = pathway.split(" ")[0]
            pathway_name = pathway.replace("{} ".format(pathway_id),"")
            fh_pathways.write("map{}\t{}\t{}\t{}\n".format(pathway_id, pathway_name, d["Pathway_category1"],
                                                           d["Pathway_category2"]))
    # Write Enzyme -> Pathway map
    with open(ec2path_out, 'w') as fh_ec2path:
        for enzyme, l in ec2path.items():
            for pathway in set(l):
                fh_ec2path.write("{}\t{}\n".format(enzyme, pathway))


def get_kegg_module_info(outdir):
    outdir = outdir.rstrip("/")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Process KEGG Module information
    logging.info("Fetching ko00002.keg from www.kegg.jp")
    url = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00002.keg&format=json"
    tmp_out = os.path.expandvars("$TMPDIR/ko00002.json")
    urllib.request.urlretrieve(url, tmp_out)
    ko2mod_out = "{}/kegg_ko2modules.tsv".format(outdir)
    modules_out = "{}/kegg_modules.tsv".format(outdir)
    with open(tmp_out) as fh:
        s = json.load(fh)
    hier = get_kegg_module_hierarchy(s)
    for f in [ko2mod_out, modules_out]:
        logging.info("Writing to {f}".format(f=f))
    with open(ko2mod_out, 'w') as fh_ko2mod, open(modules_out, 'w') as fh_modules:
        fh_modules.write("Module_id\tModule_name\tModule_category1\tModule_category2\tModule_category3\n")
        for module_name, d in hier.items():
            module_key = module_name.split(" ")[0]
            module_name = " ".join(module_name.split(" ")[1:]).lstrip()
            fh_modules.write("{}\t{}\t{}\t{}\t{}\n".format(module_key, module_name, d["Module_category1"],
                                                           d["Module_category2"], d["Module_category3"]))
            for ko in set(d["KOs"]):
                fh_ko2mod.write("{}\t{}\n".format(ko,module_key))

# Parse functions
#################


def make_orf2ko_map(df):
    orfs = []
    kos = []
    for i in df.index:
        orf = df.loc[i,"orf"]
        if df.loc[i,"ko"]==df.loc[i,"ko"]:
            for ko in df.loc[i,"ko"].split(","):
                kos.append(ko)
                orfs.append(orf)
        else:
            orfs.append(orf)
            kos.append(df.loc[i,"ko"])
    dataframe = pd.DataFrame({"orf":orfs,"ko":kos})
    return dataframe


def parse_ko_annotations(annotations, dldir, outdir):
    logging.info("Reading annotations from {}".format(annotations))
    annot = pd.read_csv(annotations, header=None, usecols=[0,2,3,6], names=["orf","evalue","score","ko"], sep="\t")
    logging.info("Loading KEGG info files from {}".format(dldir))
    ko2ec = pd.read_csv("{}/kegg_ko2ec.tsv".format(dldir), header=None, names=["ko","ec"], index_col=0, sep="\t")
    ko2path = pd.read_csv("{}/kegg_ko2pathways.tsv".format(dldir), header=None, names=["ko","pathway"], index_col=0, sep="\t")
    ko2module = pd.read_csv("{}/kegg_ko2modules.tsv".format(dldir), index_col=0, header=None, names = ["ko","module"], sep="\t")
    kos = pd.read_csv("{}/kegg_kos.tsv".format(dldir), index_col=0, sep="\t")
    modules = pd.read_csv("{}/kegg_modules.tsv".format(dldir), index_col=0, sep="\t")
    pathways = pd.read_csv("{}/kegg_pathways.tsv".format(dldir), index_col=0, sep="\t")
    orftable = make_orf2ko_map(annot)
    # Because each orf can have multiple enzyme annotations it might get placed into the same pathway several times
    # select the first combination for each orf to avoid redundancy.
    # Add KEGG Ortholog names
    orf2ko = pd.merge(orftable, kos, left_on="ko", right_index=True)
    orf2ko.set_index("orf",inplace=True)
    orf2ko.sort_index(inplace=True)
    # Map enzymes
    orf2ec = pd.merge(orftable,ko2ec,left_on="ko",right_index=True)
    # Get the first instance of each orf->ec mapping
    orf2ec = orf2ec.groupby(["orf","ec"]).first().reset_index()
    orf2ec.drop("ko", axis=1, inplace=True)
    orf2ec.set_index("orf", inplace=True)
    orf2ec.sort_index(inplace=True)
    # Map pathways
    orf2path = pd.merge(orftable,ko2path,left_on="ko",right_index=True)
    orf2path = orf2path.groupby(["orf","pathway"]).first().reset_index()
    orf2path.drop("ko", axis=1, inplace=True)
    orf2path = pd.merge(orf2path,pathways,left_on="pathway",right_index=True)
    orf2path.set_index("orf", inplace=True)
    orf2path.sort_index(inplace=True)
    # Map Modules
    orf2module = pd.merge(orftable,ko2module,left_on="ko",right_index=True)
    orf2module = orf2module.groupby(["orf","module"]).first().reset_index()
    orf2module.drop("ko", axis=1, inplace=True)
    orf2module = pd.merge(orf2module,modules,left_on="module",right_index=True)
    orf2module.set_index("orf", inplace=True)
    orf2module.sort_index(inplace=True)
    # Write files
    logging.info("Writing files to {}".format(outdir))
    orf2ko.to_csv("{}/kos.parsed.tab".format(outdir), index=True, sep="\t")
    orf2ec.to_csv("{}/enzymes.parsed.tab".format(outdir), index=True, sep="\t")
    orf2path.to_csv("{}/pathways.parsed.tab".format(outdir), index=True, sep="\t")
    orf2module.to_csv("{}/modules.parsed.tab".format(outdir), index=True, sep="\t")


# Quantify functions
####################

def process_and_sum(q_df, annot_df):
    # Merge annotations and abundance, keep ORFs without annotation as "Unclassified"
    annot_q_df = pd.merge(annot_df, q_df, left_index = True, right_index = True, how="right")
    annot_q_df.fillna("Unclassified", inplace=True)
    feature_cols = annot_df.columns
    annot_q_sum = annot_q_df.groupby(list(feature_cols)).sum().reset_index()
    annot_q_sum.set_index(feature_cols[0], inplace=True)
    return annot_q_sum


def sum_to_features(abundance, parsed):
    parsed_df = pd.read_csv(parsed, index_col=0, sep="\t")
    abundance_df = pd.read_csv(abundance, index_col=0, sep="\t")
    feature_sum = process_and_sum(abundance_df, parsed_df)
    return feature_sum


def normalize(q_df, parsed, normalize_file):
    info_df = pd.read_csv(normalize_file, header=None, sep="\t")
    info_norm_df = info_df.groupby(1).count()
    info_norm_df.columns = ["norm_factor"]
    annot_df = pd.read_csv(parsed, index_col=0, sep="\t")
    annot_cols = list(set(annot_df.columns).intersection(set(q_df.columns)))
    sample_cols = list(set(q_df.columns).difference(annot_cols))
    q_df = pd.merge(q_df,info_norm_df, left_index=True, right_index=True, how="left")
    q_df.fillna(1,inplace=True)
    q_df.loc[:,sample_cols] = q_df[sample_cols].div(q_df["norm_factor"],axis=0)
    q_df.drop("norm_factor", axis=1, inplace=True)
    return q_df


# Parser default functions
##########################


def download(args):
    logging.info("Downloading KEGG module info to {}".format(args.dldir))
    get_kegg_module_info(args.dldir)
    logging.info("Downloading KEGG ortholog info to {}".format(args.dldir))
    get_kegg_ortholog_info(args.dldir)


def parse(args):
    for f in ["kegg_ec2pathways.tsv", "kegg_ko2ec.tsv", "kegg_ko2modules.tsv", "kegg_ko2pathways.tsv",
              "kegg_kos.tsv", "kegg_modules.tsv", "kegg_pathways.tsv"]:
        if not os.path.exists("{}/{}".format(args.dldir, f)):
            download(args)
            break
    parse_ko_annotations(args.annotations, args.dldir, args.outdir)


def quant(args):
    feature_sum = sum_to_features(args.abundance, args.parsed)
    if args.normalize:
        if not os.path.exists(args.normalize):
            dldir = os.path.dirname(args.normalize)
            get_kegg_module_info(dldir)
            get_kegg_ortholog_info(dldir)
        feature_sum = normalize(feature_sum, args.parsed, args.normalize)
    feature_sum.to_csv(args.outfile, sep="\t")


def main():
    parser = ArgumentParser()
    subparser = parser.add_subparsers(title="Subcommands")
    # Download parser
    download_parser = subparser.add_parser("download", help="Download KEGG info files")
    download_parser.add_argument("dldir",
                                 help="Write files to this directory. Will be created if missing.")
    download_parser.set_defaults(func=download)
    # Annotation parser
    annot_parser = subparser.add_parser("parse", help="Parse emapper annotation files")
    annot_parser.add_argument("dldir",
                              help="Directory used to store KEGG info files (from download command)")
    annot_parser.add_argument("annotations",
                              help="emapper.py annotation files (Typically *.emapper.annotations)")
    annot_parser.add_argument("outdir",
                              help="Write parsed results to outfile")
    annot_parser.set_defaults(func=parse)
    # Sum annotations parser
    quant_parser = subparser.add_parser("quantify", help="Add abundances and summarize KEGG features")
    quant_parser.add_argument("abundance",
                              help="Abundance table")
    quant_parser.add_argument("parsed",
                              help="Parsed file (from parse command)")
    quant_parser.add_argument("outfile",
                              help="Output file with summed abundances for features")
    quant_parser.add_argument("--normalize",
                              help="Specify either kegg_ko2pathways.tsv or kegg_ko2modules.tsv to normalize abundances"
                                   "by size of the pathway/module.")
    quant_parser.set_defaults(func=quant)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
