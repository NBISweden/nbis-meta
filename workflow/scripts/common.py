#!/usr/bin/env python

import pandas as pd
import os
from os.path import join as opj, basename as bn
from snakemake.io import expand


# configuration

def prepost_string(config):
    """
    Generate preprocess string based on configuration

    :param config: Snakemake config dictionary
    :return: PREPROCESS, POSTPROCESS tuple
    """
    PREPROCESS = ""
    POSTPROCESS = ""

    preprocess_suffices = {"sortmerna": "", "trimming": "", "phixfilt": "",
                           "fastuniq": ""}

    # SortMeRNA rRNA filtering
    if config["preprocessing"]["sortmerna"]:
        PREPROCESS += ".sortmerna"
        preprocess_suffices["trimming"] = ".sortmerna"

    # Trimming
    if config["preprocessing"]["trimmomatic"]:
        PREPROCESS += ".trimmomatic"
        preprocess_suffices["phixfilt"] = preprocess_suffices[
                                              "trimming"] + ".trimmomatic"
    elif config["preprocessing"]["cutadapt"]:
        PREPROCESS += ".cutadapt"
        preprocess_suffices["phixfilt"] = preprocess_suffices[
                                              "trimming"] + ".cutadapt"
    else:
        preprocess_suffices["phixfilt"] = preprocess_suffices["trimming"]

    # Filtering
    if config["preprocessing"]["phix_filter"]:
        preprocess_suffices["fastuniq"] = preprocess_suffices[
                                              "phixfilt"] + ".phixfilt"
        PREPROCESS += ".phixfilt"
    else:
        preprocess_suffices["fastuniq"] = preprocess_suffices["phixfilt"]

    # Deduplication
    if config["preprocessing"]["fastuniq"]:
        PREPROCESS += ".fastuniq"

    if PREPROCESS != "":
        config["preprocess"] = True
    else:
        config["preprocess"] = False

    if config["remove_duplicates"]:
        POSTPROCESS += ".markdup"

    return PREPROCESS, POSTPROCESS, preprocess_suffices, config


def parse_samples(df, config, PREPROCESS):
    assemblies = {}
    samples = {}
    df.fillna("", inplace=True)

    for i in list(df.index):
        # Add sample to dict
        sample = df.iloc[i]["sample"]
        if sample not in samples.keys():
            samples[sample] = {}
        # Add unit to dict
        unit = str(df.iloc[i]["unit"])
        if unit not in samples[sample].keys():
            samples[sample][unit] = {}
        R1 = df.iloc[i]["fq1"]
        groups = []
        r2 = False

        # Set preprocessed file paths
        R1_p = opj(config["paths"]["results"], "intermediate", "preprocess",
                   "{}_{}_R1{}.fastq.gz".format(sample, unit, PREPROCESS))
        R2_p = opj(config["paths"]["results"], "intermediate", "preprocess",
                   "{}_{}_R2{}.fastq.gz".format(sample, unit, PREPROCESS))
        se_p = opj(config["paths"]["results"], "intermediate", "preprocess",
                   "{}_{}_se{}.fastq.gz".format(sample, unit, PREPROCESS))

        # Initiate keys for all assembly group values
        if "assembly" in df.columns:
            assem_list = df.iloc[i]["assembly"].split(",")
            for a in assem_list:
                if a not in assemblies.keys() and a != "":
                    assemblies[a] = {}
        # Handling of paired and/or single end sequence files If the sample
        # annotation file has a 'pair' column, add the read files as 'R1' and
        # 'R2'
        if "fq2" in df.columns and df.iloc[i]["fq2"]:
            R2 = df.iloc[i]["fq2"]
            r2 = True
            samples[sample][unit]["R1"] = R1
            samples[sample][unit]["R2"] = R2
            # Add filepaths to preprocessed output files for each of the read
            # files in each of the assemblies. This will be the initial
            # input to the assembly rule
            for a in assemblies:
                if sample not in assemblies[a].keys():
                    assemblies[a][sample] = {unit: {}}
                if r2:
                    assemblies[a][sample][unit]["R1"] = [R1_p]
                    assemblies[a][sample][unit]["R2"] = [R2_p]
                else:
                    assemblies[a][sample][unit]["se"] = [se_p]
        # If there is no 'fq2' column, add the single file path as 'se'
        else:
            samples[sample][unit]["se"] = R1
            for a in assemblies:
                if sample not in assemblies[a].keys():
                    assemblies[a][sample] = {unit: {}}
                assemblies[a][sample][unit]["se"] = [se_p]
    return samples, assemblies


def check_uppmax(config):
    """
    Set specific params for running on uppmax

    :param config: Snakemake config
    :return: updated config dictionary
    """
    import platform
    hostname = platform.node()
    if 'uppmax.uu.se' in hostname:
        config["runOnUppMax"] = True
        # Set temp to $TMPDIR
        config["paths"]["temp"] = "$TMPDIR"
    else:
        config["runOnUppMax"] = False
    return config


def check_annotation(config):
    """
    Checks whether to run annotation/assembly

    :param config: Snakemake config
    :return: Updated config dict
    """
    # Check whether to set annotation downstream of assembly
    tools = [config["pfam"], config["taxonomic_annotation"], config["infernal"],
             config["eggnog"], config["rgi"]]
    if True in tools:
        config["run_annotation"] = True
        # if True also assume the user wants assembly
        config["run_assembly"] = True
        # Set megahit as default unless metaspades is set
        if not config["assembly"]["megahit"] and not config["assembly"]["metaspades"]:
            config["assembly"]["megahit"] = True
    else:
        config["annotation"] = False
    return config


def check_assembly(config, assemblies):
    """
    Check assemblies and config settings

    :param config: Snakemake config
    :param assemblies: Assembly dictionary
    :return: Tuple of updated config and assembly dict
    """
    if len(assemblies) > 0:
        if config["metaspades"]:
            # Remove single-end only assemblies
            # that Metaspades won't be able to run
            assemblies = filter_metaspades_assemblies(assemblies)
        # Add information on binning
        binning = False
        if config["maxbin"] or config["concoct"] or config["metabat"]:
            binning = True
        if not binning:
            config["checkm"] = False
            config["gtdbtk"] = False
    return config, assemblies


def check_classifiers(config):
    """
    Set paths and params specific to classifiers

    :param config: Snakemake config
    :return: Updated config dict
    """
    # Add read-based config info
    config["centrifuge_index_path"] = config["centrifuge_base"] = config[
        "centrifuge_dir"] = ""
    if config["centrifuge"]:
        # Check if custom database exists
        custom = expand("{b}.{i}.cf", b=config["centrifuge_custom"],
                        i=[1, 2, 3])
        if list(set([os.path.exists(x) for x in custom]))[0]:
            config["centrifuge_index_path"] = config["centrifuge_custom"]
        # If not, use prebuilt default
        else:
            config["centrifuge_index_path"] = "resources/centrifuge/{}".format(
                config["centrifuge_prebuilt"])
        # Set centrifuge index config variables
        config['centrifuge_dir'] = os.path.dirname(
            config['centrifuge_index_path'])
        config['centrifuge_base'] = bn(config['centrifuge_index_path'])

    config["kraken_index_path"] = config["kraken_params"] = ""
    if config["kraken"]:
        # Check if custom database exists
        custom = expand(opj(config["kraken_custom"], "{n}.k2d"),
                        n=["hash", "opts", "taxo"])
        if list(set(os.path.exists(x) for x in custom))[0]:
            config["kraken_index_path"] = config["kraken_custom"]
        # If not, use prebuilt or standard
        elif config["kraken_standard_db"]:
            config["kraken_index_path"] = opj("resources", "kraken",
                                              "standard")
        else:
            config["kraken_index_path"] = opj("resources", "kraken",
                                              "prebuilt",
                                              config["kraken_prebuilt"])
        if config["kraken_reduce_memory"]:
            config["kraken_params"] = "--memory-mapping"
    return config


# preprocessing functions

def preprocessing_input(config):
    if config["preprocess"] or config["preprocessing"]["fastqc"]:
        return opj(config["paths"]["results"], "report", "samples_report.html")
    return []


def link(target, link_name):
    """
    Generates symlinks with absolute paths

    :param target:
    :param link_name:
    :return:
    """
    src = os.path.abspath(target)
    dst = os.path.abspath(link_name)
    os.symlink(src, dst)


def is_pe(d):
    if "R2" in d.keys():
        return True
    return False


def get_all_files(samples, dir, suffix="", nested=False):
    """
    Returns a list of files based on samples and directory

    :param samples: Samples dictionary
    :param dir: Directory to find files in
    :param suffix: Suffix of files to return
    :param nested: If True look for files inside sample_run directory
    :return:
    """
    files = []
    for sample in samples:
        for unit in samples[sample].keys():
            if nested:
                d = "{}/{}_{}".format(dir, sample, unit)
            else:
                d = "{}".format(dir)
            if is_pe(samples[sample][unit]):
                files.append(opj(d, "{}_{}_pe{}".format(sample, unit, suffix)))
            else:
                files.append(opj(d, "{}_{}_se{}".format(sample, unit, suffix)))
    return files


def multiqc_input(samples, config):
    files = []
    pre, post, d, _ = prepost_string(config)
    for sample in samples.keys():
        for unit in samples[sample].keys():
            if is_pe(samples[sample][unit]):
                pairs = ["R1", "R2"]
                seq_type = "pe"
            else:
                pairs = ["se"]
                seq_type = "se"
            files += get_fastqc_files(sample, unit, pairs, config, pre)
            files += get_trim_logs(sample, unit, pairs, config, d)
            files += get_filt_logs(sample, unit, seq_type, config, d)
            files += get_sortmerna_logs(sample, unit, seq_type, config)
    return files


def get_fastqc_files(sample, unit, pairs, config, pre):
    """Get all fastqc output"""
    if config["preprocessing"]["fastqc"]:
        files = expand(opj(config["paths"]["results"], "intermediate", "fastqc",
                           "{sample}_{unit}_{pair}{PREPROCESS}_fastqc.zip"),
                       sample=sample, unit=unit, pair=pairs, PREPROCESS=pre)
        return files
    return []


def get_trim_logs(sample, unit, pairs, config, d):
    if not config["preprocessing"]["trimmomatic"] and not config["preprocessing"]["cutadapt"]:
        return []
    if config["preprocessing"]["trimmomatic"]:
        trimmer = "trimmomatic"
    else:
        trimmer = "cutadapt"
    files = expand(opj(config["paths"]["results"], "intermediate", "preprocess",
                       "{sample}_{unit}_{pair}{s}.{trimmer}.log"),
                   sample=sample, unit=unit, pair=pairs, s=d["trimming"],
                   trimmer=trimmer)
    return files


def get_filt_logs(sample, unit, seq_type, config, d):
    if not config["preprocessing"]["phix_filter"]:
        return []
    files = expand(opj(config["paths"]["results"], "intermediate", "preprocess",
                       "{sample}_{unit}_PHIX_{seq_type}{s}.log"), sample=sample,
                   unit=unit, seq_type=seq_type, s=d["phixfilt"])
    return files


def get_sortmerna_logs(sample, unit, seq_type, config):
    if not config["preprocessing"]["sortmerna"]:
        return []
    files = expand(opj(config["paths"]["results"], "intermediate", "preprocess",
                       "{sample}_{unit}_{seq_type}.sortmerna.log"),
                   sample=sample, unit=unit, seq_type=seq_type)
    return files


def get_trimmomatic_string(seq_type, config):
    """
    Generates trimsetting string for Trimmomatic

    :param seq_type: PE or SE depending on sequencing type
    :return: string
    """
    trim_adapters = config["trimmomatic"]["trim_adapters"]
    adapter_fasta_dir = "$CONDA_PREFIX/share/trimmomatic/adapters"
    # Get params based on sequencing type
    param_dict = config["trimmomatic"][seq_type]
    # Set path to adapter
    adapter = "{}/{}.fa".format(adapter_fasta_dir, param_dict["adapter"])
    adapter_params = param_dict["adapter_params"]
    pre_adapter_params = param_dict["pre_adapter_params"]
    post_adapter_params = param_dict["post_adapter_params"]
    trimsettings = pre_adapter_params
    if trim_adapters:
        trimsettings = " {} ILLUMINACLIP:{}:{}".format(pre_adapter_params,
                                                       adapter, adapter_params)
    return "{} {}".format(trimsettings, post_adapter_params)


def get_sortmerna_ref_string(dbs):
    """
    Constructs the SortMeRNA --ref string

    :param dbs: Sortmerna databases from config
    :return: STRING, STRING formatted string
    """
    files = ["resources/rRNA_databases/{db}".format(db=db) for db in dbs]
    ref_string = ":".join(["{},{}".format(f, f) for f in files])
    return ref_string


# assembly functions

def filter_metaspades_assemblies(d):
    """
    This function removes assemblies that contain only single-end samples

    :param d: Assembly group dictionary
    :return: Dictionary containing only assemblies with at least 1 paired sample
    """
    se_only = []
    for assembly in d.keys():
        i = 0
        for sample in d[assembly].keys():
            for unit in d[assembly][sample].keys():
                if is_pe(d[assembly][sample][unit]):
                    i += 1
                    break
        if i == 0:
            se_only.append(assembly)
    for assembly in se_only:
        del d[assembly]
    return d


def get_all_group_files(assembly_dict):
    files = []
    for sample in assembly_dict.keys():
        for unit in assembly_dict[sample].keys():
            for pair in assembly_dict[sample][unit].keys():
                files.append(assembly_dict[sample][unit][pair][0])
    return files


def get_bamfiles(g, assembly_dict, results_path, POSTPROCESS):
    files = []
    for sample in assembly_dict.keys():
        for unit in assembly_dict[sample].keys():
            if "R2" in assembly_dict[sample][unit].keys():
                seq_type = "pe"
            else:
                seq_type = "se"
            files.append(opj(results_path, "assembly", g, "mapping",
                             "{}_{}_{}{}.bam".format(sample, unit, seq_type,
                                                     POSTPROCESS)))
    return files


def rename_records(f, fh, i):
    """
    Prepends a number to read ids

    :param f: Input fastq file (gzipped)
    :param fh: Output filehandle
    :param i: File index to prepend to reads
    :return: Output filehandle
    """
    from Bio import SeqIO
    import gzip as gz
    for record in SeqIO.parse(gz.open(f, 'rt'), 'fastq'):
        record.id = "{}_{}".format(i, record.id)
        SeqIO.write(record, fh, "fastq")
    return fh


# binning functions

def binning_input(config, assemblies):
    """
    Generates input list for the binning part of the workflow

    :param config: Snakemake config
    :param assemblies: Dictionary of assemblies
    :return:
    """
    binners = get_binners(config)
    bin_input = expand(
        opj(config["paths"]["results"], "binning", "{binner}", "{group}", "{l}",
            "summary_stats.tsv"), binner=binners, group=assemblies.keys(),
        l=config["min_contig_length"])

    if config["checkm"]:
        bin_input.append(opj(config["paths"]["results"], "report", "checkm", "checkm.stats.tsv"))
        bin_input.append(
            opj(config["paths"]["results"], "report", "checkm", "checkm.profiles.tsv"))
    if config["gtdbtk"]:
        bin_input.append(opj(config["paths"]["results"], "report", "gtdbtk", "gtdbtk.summary.tsv"))
        bin_input.append(
            opj(config["paths"]["results"], "report", "bin_annotation", "tRNA.total.tsv"))
        bin_input.append(
            opj(config["paths"]["results"], "report", "bin_annotation", "rRNA.types.tsv"))
    return bin_input


def get_fw_reads(config, samples, p):
    """
    MaxBin2 only uses unpaired reads for mapping with bowtie2.
    Here we iterate over all samples
    """
    files = []
    for sample in samples.keys():
        for unit in samples[sample].keys():
            if "R1" in samples[sample][unit].keys():
                f = opj(config["paths"]["results"], "intermediate", "preprocess",
                        "{sample}_{unit}_R1{p}.fastq.gz".format(sample=sample,
                                                                unit=unit, p=p))
            else:
                f = opj(config["paths"]["results"], "intermediate", "preprocess",
                        "{sample}_{unit}_se{p}.fastq.gz".format(sample=sample,
                                                                unit=unit, p=p))
            files.append(f)
    reads_string = ""
    for i, f in enumerate(files, start=1):
        reads_string += "-reads{i} {f} ".format(i=i, f=f)
    return reads_string


def get_tree_settings(config):
    """
    Return checkm parameter based on tree settings

    :param config:
    :return:
    """
    if config["checkm_reduced_tree"]:
        return "-r"
    return ""


def get_binners(config):
    """
    Return a list of binners used
    :param config:
    :return:
    """
    binners = []
    if config["metabat"]:
        binners.append("metabat")
    if config["concoct"]:
        binners.append("concoct")
    if config["maxbin"]:
        binners.append("maxbin")
    return binners


def get_fields(f):
    """
    Extract assembly, binner and length fields from file path
    :param f: Input file
    :return:
    """
    items = f.split("/")
    return items[-3], items[-4], items[-5]


def assign_fields(x, l, group, binner):
    """
    Assign assembly, binner and length fields
    :param x: pandas DataFrame
    :param l: minimum contig length used
    :param group: assembly group
    :param binner: binner used
    :return: updated pandas DataFrame
    """
    rows = x.shape[0]
    x = x.assign(binner=pd.Series([binner] * rows, index=x.index))
    x = x.assign(min_contig_length=pd.Series([l] * rows, index=x.index))
    x = x.assign(assembly=pd.Series([group] * rows, index=x.index))
    return x


def concatenate(input):
    """
    Concatenate bin info dataframes
    :param input:
    :return:
    """
    df = pd.DataFrame()
    for i, f in enumerate(input):
        l, group, binner = get_fields(f)
        _df = pd.read_csv(f, sep="\t", index_col=0)
        rows = _df.shape[0]
        if rows == 0:
            continue
        _df = assign_fields(_df, l, group, binner)
        if i == 0:
            df = _df.copy()
        else:
            _df = _df.loc[:, df.columns]
            df = pd.concat([df, _df])
    return df


# annotation functions

def annotation_input(config, assemblies):
    input = []
    for group in assemblies.keys():
        # Add orfcalling results
        input.append(opj(config["paths"]["results"], "annotation", group,
                         "final_contigs.gff"))
        if config["infernal"]:
            input.append(opj(config["paths"]["results"], "annotation", group,
                             "final_contigs.cmscan"))
        if config["tRNAscan"]:
            input.append(
                opj(config["paths"]["results"], "annotation", group, "tRNA.out"))
        # Add EGGNOG annotation
        if config["eggnog"]:
            input += expand(opj(config["paths"]["results"], "annotation", group,
                                "{db}.parsed.{fc}.tsv"),
                            db=["enzymes", "pathways", "kos", "modules"],
                            fc=["raw", "tpm"])
        # Add PFAM annotation
        if config["pfam"]:
            input += expand(opj(config["paths"]["results"], "annotation", group,
                                "pfam.parsed.{fc}.tsv"), fc=["tpm", "raw"])
        # Add taxonomic annotation
        if config["taxonomic_annotation"]:
            input += expand(
                opj(config["paths"]["results"], "annotation", group, "taxonomy",
                    "tax.{fc}.tsv"), fc=["tpm", "raw"])
        # Add Resistance Gene Identifier output
        if config["rgi"]:
            input += expand(opj(config["paths"]["results"], "annotation", group,
                                "rgi.{fc}.tsv"), fc=["raw", "tpm"])
            input.append(
                opj(config["paths"]["results"], "annotation", group, "rgi.out.txt"))
    return input


def parse_cmout(f):
    with open(f) as fh:
        lines = []
        idnums = {}
        for i, line in enumerate(fh):
            if line[0] == "#": continue
            line = line.rstrip()
            line = line.rsplit()
            indices = [1, 2, 3, 5, 9, 10, 11, 14, 16, 17]
            target_name, target_accession, query, clan, start, end, strand, gc, score, evalue = [
                line[x] for x in indices]
            try:
                idnum = idnums[query]
            except KeyError:
                idnum = 0
                idnums[query] = idnum
            this_idnum = idnum + 1
            idnums[query] = this_idnum
            attributes = ["ID=" + query + "ncRNA_" + str(this_idnum),
                          "Name=" + target_name,
                          "Accession=" + target_accession, "Clan=" + clan,
                          "GC=" + gc, "E-value=" + evalue]
            # seqid, source, type, start, end, score, strand, phase, attributes
            gffline = " ".join(
                [query, "cmscan", "ncRNA", start, end, score, strand, ".",
                 ";".join(attributes)])
            lines.append(gffline)
    return lines


## Quantification functions

def markdup_mem(wildcards, cores):
    """
    Calculates the memory to allocate when running MarkDuplicates
    :param wildcards:
    :param cores: number of cores for currently running workflow
    :return:
    """
    threads = min(cores, 10)
    mem_gb_per_thread = 2
    return int(mem_gb_per_thread*threads)


def get_fc_files(wildcards, file_type):
    g = wildcards.group
    files = []
    for sample in assemblies[g].keys():
        for unit in assemblies[g][sample].keys():
            if "se" in assemblies[g][sample][unit].keys():
                files.append(
                    opj(config["paths"]["results"], "assembly", g, "mapping",
                        sample + "_" + unit + "_se.fc.{}.tsv".format(
                            file_type)))
            else:
                files.append(
                    opj(config["paths"]["results"], "assembly", g, "mapping",
                        sample + "_" + unit + "_pe.fc.{}.tsv".format(
                            file_type)))
    return files


# classification functions

def classify_input(config):
    f = []
    if config["kraken"]:
        f.append(opj(config["paths"]["results"], "report", "kraken", "kraken.krona.html"))
    if config["metaphlan"]:
        f.append(opj(config["paths"]["results"], "report", "metaphlan", "metaphlan.html"))
    if config["centrifuge"]:
        f.append(
            opj(config["paths"]["results"], "report", "centrifuge", "centrifuge.krona.html"))
    return f


def krona_input(config, samples, classifier):
    input_string = ""
    files = get_all_files(samples, opj(config["paths"]["results"], classifier),
                          ".kreport")
    for f in files:
        sample_unit = bn(f).replace("_pe.kreport", "").replace("_se.kreport",
            "")
        input_string += " {},{}".format(f, sample_unit)
    return input_string


def get_centrifuge_index_url(config):
    url_base = "ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data"
    d = {'nt_2018_2_12': 'nt_2018_2_12.tar.gz',
         'nt_2018_3_3': 'nt_2018_3_3.tar.gz', 'p+h+v': 'p+h+v.tar.gz',
         'p_compressed+h+v': 'p_compressed+h+v.tar.gz',
         'p_compressed_2018_4_15': 'p_compressed_2018_4_15.tar.gz'}
    try:
        url = "{}/{}".format(url_base, d[config["centrifuge_base"]])
    except KeyError:
        url = ""
    return url


def get_kraken_index_url(kraken_prebuilt, version=False):
    """
    Downloads latest prebuilt kraken index

    :param kraken_prebuilt: prebuilt version
    :param version: Return db version on True
    :return: url text string
    """
    import subprocess
    import re
    url_base = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs"
    r = subprocess.run(["curl", "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/"],
                       capture_output=True)
    _files = [x.split(" ")[-1] for x in r.stdout.decode().split("\n")]
    # List all tar archives
    files = [x for x in _files if ".tgz" in x]
    # Figure out versions
    versions = []
    types = []
    for t in ["_".join(x.split("_")[0:2]) for x in files]:
        if t == "minikraken_8GB":
            types.append(t)
            versions.append("")
            continue
        v = re.search("\d+.*", t.split("_")[-1]).group()
        types.append(t.replace(v, ""))
        versions.append(v)
    # Build a DataFrame and sort by version
    df = pd.DataFrame({'type': types, 'version': versions, 'file': files})
    df = df.sort_values("version", ascending=False)
    if df.shape[0] == 0:
        return ""
    f = df.loc[df.type == kraken_prebuilt].head(1)["file"].values[0]
    if version:
        return os.path.splitext(f)[0]
    return "{}/{}".format(url_base, f)


def metaphlan_krona_string(input):
    """
    Returns the input string (<file1>,<name1> ... <fileN>,<nameN>) for krona
    :param input: krona-ready metaphlan tables
    :return: krona-formatted input string
    """
    s = []
    for f in input:
        name = bn(f).replace("_pe.krona", "").replace("_se.krona", "")
        s.append("{},{}".format(f, name))
    return " ".join(s)
