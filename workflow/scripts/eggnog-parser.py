#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import json
import os
from urllib import request
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
                    #kos = module['children']
                    hier[module_name] = {"Module_category1": c1, "Module_category2": c2,
                                         "Module_category3": c3}
                    #for ko in kos:
                    #    ko_name = ko['name']
                    #    ko = ko_name.split(" ")[0]
                    #    hier[module_name]["KOs"].append(ko)
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


def setup_tmpdir(tmpdir):
    if not os.path.isdir(os.path.expandvars(tmpdir)):
        try:
            os.makedirs(tmpdir, exist_ok=True)
        except PermissionError:
            tmpdir = "temp"
            os.makedirs(tmpdir, exist_ok=True)


def get_kegg_ortholog_info(outdir, tmpdir="$TMPDIR"):
    outdir = outdir.rstrip("/")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    setup_tmpdir(tmpdir)
    url = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json"
    logging.info("Fetching ko00001.keg from www.kegg.jp")
    # Download file
    tmp_out = os.path.expandvars(f"{tmpdir}/ko00001.json")
    request.urlretrieve(url, tmp_out)
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
            pathway_name = pathway.replace("{} ".format(pathway_id), "")
            fh_pathways.write("map{}\t{}\t{}\t{}\n".format(pathway_id, pathway_name, d["Pathway_category1"],
                                                           d["Pathway_category2"]))
    # Write Enzyme -> Pathway map
    with open(ec2path_out, 'w') as fh_ec2path:
        for enzyme, l in ec2path.items():
            for pathway in set(l):
                fh_ec2path.write("{}\t{}\n".format(enzyme, pathway))


def get_kegg_module_info(outdir, tmpdir="$TMPDIR"):
    outdir = outdir.rstrip("/")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    setup_tmpdir(tmpdir)
    outdir = outdir.rstrip("/")
    # Process KEGG Module information
    logging.info("Fetching ko00002.keg from www.kegg.jp")
    url = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00002.keg&format=json"
    tmp_out = os.path.expandvars(f"{tmpdir}/ko00002.json")
    request.urlretrieve(url, tmp_out)
    modules_out = "{}/kegg_modules.tsv".format(outdir)
    with open(tmp_out) as fh:
        s = json.load(fh)
    hier = get_kegg_module_hierarchy(s)
    logging.info("Writing to {}".format(modules_out))
    with open(modules_out, 'w') as fh_modules:
        fh_modules.write("Module_id\tModule_name\tModule_category1\tModule_category2\tModule_category3\n")
        for module_name, d in hier.items():
            module_key = module_name.split(" ")[0]
            module_name = " ".join(module_name.split(" ")[1:]).lstrip()
            fh_modules.write("{}\t{}\t{}\t{}\t{}\n".format(module_key, module_name, d["Module_category1"],
                                                           d["Module_category2"], d["Module_category3"]))


# Parse functions
#################
def feature2orf(df, feature, lstrip=False, match=False, extra=""):
    orfs = []
    features = []
    extras = []
    for i in df.index:
        orf = df.loc[i,"orf"]
        feats = df.loc[i,feature].split(",")
        if lstrip:
            feats = [x.lstrip(lstrip) for x in feats]
        if match:
            feats = [x for x in feats if match in x]
        if extra != "":
            extras.append(df.loc[i, extra])
        orfs+=[orf]*len(feats)
        features+=feats
    if extra == "":
        return pd.DataFrame([orfs, features], index=["orf", feature]).T
    else:
        return pd.DataFrame([orfs, features, extras], index=["orf", feature, extra]).T

def parse_ko_annotations(annotations, dldir, outdir, map_go=False):
    logging.info("Reading annotations from {}".format(annotations))
    annot = pd.read_csv(annotations, header=None, usecols=[0, 2, 3, 6, 8, 9, 10, 14, 15, 21],
                        names=["orf", "evalue", "score", "GO", "ko", "pathway", "module", "tc", "cazy", "desc"],
                        sep="\t")
    logging.info("Mapping ORFs to KEGG Orthologs")
    orf2ko = feature2orf(annot.loc[annot["ko"] == annot["ko"]], "ko", lstrip="ko:")
    logging.info("{} ORFs with {} KOs".format(len(orf2ko["orf"].unique()), len(orf2ko["ko"].unique())))
    logging.info("Mapping ORFs to modules")
    orf2module = feature2orf(annot.loc[annot["module"]==annot["module"]], "module")
    logging.info("{} ORFs in {} modules".format(len(orf2module["orf"].unique()), len(orf2module["module"].unique())))
    logging.info("Mapping ORFs to pathways")
    orf2pathway = feature2orf(annot.loc[annot["pathway"]==annot["pathway"]], "pathway", match="map")
    logging.info("{} ORFs in {} Pathways".format(len(orf2pathway["orf"].unique()), len(orf2pathway["pathway"].unique())))
    logging.info("Mapping ORFs to transporters")
    orf2tc = feature2orf(annot.loc[annot["tc"]==annot["tc"]], feature="tc")
    logging.info("{} ORFs in {} transporters".format(len(orf2tc["orf"].unique()), len(orf2tc["tc"].unique())))
    orf2tc.set_index("orf", inplace=True)
    orf2tc.sort_index(inplace=True)
    logging.info("Mapping ORFs to CAZY")
    orf2caz = feature2orf(annot.loc[annot["cazy"]==annot["cazy"]], feature="cazy")
    logging.info("{} ORFs in {} CAZYs".format(len(orf2caz["orf"].unique()), len(orf2caz["cazy"].unique())))
    orf2caz.set_index("orf", inplace=True)
    orf2caz.sort_index(inplace=True)

    logging.info("Loading KEGG info files from {}".format(dldir))
    ko2ec = pd.read_table("{}/kegg_ko2ec.tsv".format(dldir), header=None, names=["ko", "ec"], index_col=0)
    kos = pd.read_table("{}/kegg_kos.tsv".format(dldir), index_col=0)
    modules = pd.read_table("{}/kegg_modules.tsv".format(dldir), index_col=0)
    pathways = pd.read_table("{}/kegg_pathways.tsv".format(dldir), index_col=0)
    # Because each orf can have multiple enzyme annotations it might get placed into the same pathway several times
    # select the first combination for each orf to avoid redundancy.
    # Add KEGG Ortholog names
    orf2ko = pd.merge(orf2ko, kos, left_on="ko", right_index=True, how="left")
    orf2ko.set_index("orf", inplace=True)
    orf2ko.sort_index(inplace=True)
    # Map enzymes
    logging.info("Mapping enzymes")
    orf2ec = pd.merge(orf2ko, ko2ec, left_on="ko", right_index=True)
    # Get the first instance of each orf->ec mapping
    orf2ec = orf2ec.groupby(["orf", "ec"]).first().reset_index()
    orf2ec.drop("ko", axis=1, inplace=True)
    orf2ec.set_index("orf", inplace=True)
    orf2ec.sort_index(inplace=True)
    # Add pathway info
    orf2pathway = pd.merge(orf2pathway, pathways, left_on="pathway", right_index=True, how="left")
    orf2pathway.set_index("orf", inplace=True)
    orf2pathway.sort_index(inplace=True)
    # Add module info
    orf2module = pd.merge(orf2module, modules, left_on="module", right_index=True, how="left")
    orf2module.set_index("orf", inplace=True)
    orf2module.sort_index(inplace=True)
    # Map GO IDs
    if map_go:
        logging.info("Mapping ORFs to GO terms")
        orf2go = feature2orf(annot.loc[annot["GO"]==annot["GO"]], feature="GO")
        logging.info("{} ORFs in {} GO terms".format(len(orf2go["orf"].unique()), len(orf2go["GO"].unique())))
    # Write files
    logging.info("Writing files to {}".format(outdir))
    orf2ko.to_csv("{}/kos.parsed.tsv".format(outdir), index=True, sep="\t")
    orf2ec.to_csv("{}/enzymes.parsed.tsv".format(outdir), index=True, sep="\t")
    orf2pathway.to_csv("{}/pathways.parsed.tsv".format(outdir), index=True, sep="\t")
    orf2module.to_csv("{}/modules.parsed.tsv".format(outdir), index=True, sep="\t")
    orf2tc.to_csv("{}/tc.parsed.tsv".format(outdir), index=True, sep="\t")
    orf2caz.to_csv("{}/cazy.parsed.tsv".format(outdir), index=True, sep="\t")
    if map_go:
        orf2go.to_csv("{}/gos.parsed.tsv".format(outdir), index=True, sep="\t")


# Quantify functions
####################

def process_and_sum(q_df, annot_df):
    # Merge annotations and abundance, keep ORFs without annotation as "Unclassified"
    annot_q_df = pd.merge(annot_df, q_df, left_index=True, right_index=True, how="right")
    annot_q_df.fillna("Unclassified", inplace=True)
    feature_cols = annot_df.columns
    annot_q_sum = annot_q_df.groupby(list(feature_cols)).sum().reset_index()
    annot_q_sum.set_index(feature_cols[0], inplace=True)
    return annot_q_sum


def sum_to_features(abundance, parsed):
    parsed_df = pd.read_table(parsed, index_col=0)
    abundance_df = pd.read_table(abundance, index_col=0)
    feature_sum = process_and_sum(abundance_df, parsed_df)
    return feature_sum


def normalize(q_df, parsed, normalize_file):
    info_df = pd.read_table(normalize_file, header=None)
    info_norm_df = info_df.groupby(1).count()
    info_norm_df.columns = ["norm_factor"]
    annot_df = pd.read_table(parsed, index_col=0)
    annot_cols = list(set(annot_df.columns).intersection(set(q_df.columns)))
    sample_cols = list(set(q_df.columns).difference(annot_cols))
    q_df = pd.merge(q_df, info_norm_df, left_index=True, right_index=True, how="left")
    q_df.fillna(1, inplace=True)
    q_df.loc[:, sample_cols] = q_df[sample_cols].div(q_df["norm_factor"], axis=0)
    q_df.drop("norm_factor", axis=1, inplace=True)
    return q_df


def merge_files(files, sum_abundance=False):
    df = pd.DataFrame()
    for i, f in enumerate(files):
        _df = pd.read_table(f, header=0, sep="\t")
        if sum_abundance:
            group_cols = list(_df.columns)[1:-1]
            if len(group_cols) == 0:
                group_cols = [list(_df.columns)[0]]
        else:
            group_cols = list(_df.columns)[0:-1]
        # Sum dataframe
        _df = _df.groupby(group_cols).sum().reset_index()
        if i == 0:
            df = _df.copy()
        else:
            df = pd.merge(df, _df, on=group_cols, how="outer")
    df.set_index(group_cols[0], inplace=True)
    df.fillna(0, inplace=True)
    return df


# Parser default functions
##########################


def download(args):
    if not os.path.exists(args.dldir):
        os.makedirs(args.dldir)
    logging.info("Downloading KEGG module info to {}".format(args.dldir))
    get_kegg_module_info(args.dldir)
    logging.info("Downloading KEGG ortholog info to {}".format(args.dldir))
    get_kegg_ortholog_info(args.dldir)


def parse(args):
    for f in ["kegg_ec2pathways.tsv", "kegg_ko2ec.tsv", "kegg_ko2pathways.tsv",
              "kegg_kos.tsv", "kegg_modules.tsv", "kegg_pathways.tsv"]:
        if not os.path.exists("{}/{}".format(args.dldir, f)):
            download(args)
            break
    parse_ko_annotations(args.annotations, args.dldir, args.outdir, args.map_go)


def quant(args):
    feature_sum = sum_to_features(args.abundance, args.parsed)
    if args.normalize:
        if not os.path.exists(args.normalize):
            dldir = os.path.dirname(args.normalize)
            get_kegg_module_info(dldir)
            get_kegg_ortholog_info(dldir)
        feature_sum = normalize(feature_sum, args.parsed, args.normalize)
    feature_sum.to_csv(args.outfile, sep="\t")


def merge(args):
    df = merge_files(files=args.files, sum_abundance=args.sum)
    df.to_csv(args.outfile, sep="\t", index=True, header=True)


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
                              help="Output directory for parsed files")
    annot_parser.add_argument("--map_go", action="store_true",
                              help="Also map GO terms (can take a long time)")
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
    # Merge files parser
    merge_parser = subparser.add_parser("merge", help="Merge two or more files from quantify step")
    merge_parser.add_argument("files", nargs='+',
                              help="Files from quantify step")
    merge_parser.add_argument("outfile", type=str,
                              help="Output file")
    merge_parser.add_argument("--sum", action="store_true",
                              help="When merging, sum abundances to second column. (default = False)")
    merge_parser.set_defaults(func=merge)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
