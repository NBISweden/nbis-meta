#!/usr/bin/env python

import pandas as pd
import sys
import os
from os.path import basename as bn
from snakemake.io import expand


# configuration


def mem_allowed(wildcards, threads):
    mem_per_core = config["mem_per_core"]
    return max(threads * mem_per_core, mem_per_core)


def prepost_string(config):
    """
    Generate preprocess string based on configuration

    :param config: Snakemake config dictionary
    :return: PREPROCESS, POSTPROCESS tuple
    """
    PREPROCESS = ""
    POSTPROCESS = ""

    preprocess_suffices = {
        "sortmerna": "",
        "trimming": "",
        "phixfilt": "",
        "fastuniq": "",
    }

    # SortMeRNA rRNA filtering
    if config["preprocessing"]["sortmerna"]:
        PREPROCESS += ".sortmerna"
        preprocess_suffices["trimming"] = ".sortmerna"

    # Trimming
    if config["preprocessing"]["trimmomatic"]:
        PREPROCESS += ".trimmomatic"
        preprocess_suffices["phixfilt"] = (
            preprocess_suffices["trimming"] + ".trimmomatic"
        )
    elif config["preprocessing"]["cutadapt"]:
        PREPROCESS += ".cutadapt"
        preprocess_suffices["phixfilt"] = preprocess_suffices["trimming"] + ".cutadapt"
    else:
        preprocess_suffices["phixfilt"] = preprocess_suffices["trimming"]

    # Filtering
    if config["preprocessing"]["phix_filter"]:
        preprocess_suffices["fastuniq"] = preprocess_suffices["phixfilt"] + ".phixfilt"
        PREPROCESS += ".phixfilt"
    else:
        preprocess_suffices["fastuniq"] = preprocess_suffices["phixfilt"]

    # Deduplication
    if config["preprocessing"]["fastuniq"]:
        PREPROCESS += ".fastuniq"

    if PREPROCESS != "":
        config["run_preprocessing"] = True
    else:
        config["run_preprocessing"] = False

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
        if "unit" not in df.columns:
            unit = "1"
        else:
            unit = str(df.iloc[i]["unit"])
        if unit not in samples[sample].keys():
            samples[sample][unit] = {}
        R1 = df.iloc[i]["fq1"]
        r2 = False

        # Set preprocessed file paths
        R1_p = config["paths"][
            "results"
        ] + "/intermediate/preprocess/{}_{}_R1{}.fastq.gz".format(
            sample, unit, PREPROCESS
        )
        R2_p = config["paths"][
            "results"
        ] + "/intermediate/preprocess/{}_{}_R2{}.fastq.gz".format(
            sample, unit, PREPROCESS
        )
        se_p = config["paths"][
            "results"
        ] + "/intermediate/preprocess/{}_{}_se{}.fastq.gz".format(
            sample, unit, PREPROCESS
        )

        # Initiate keys for all assembly group values
        if "assembly" in df.columns:
            assem_list = df.iloc[i]["assembly"].split(",")
            assem_list = [a for a in assem_list if a != ""]
            for a in assem_list:
                if a not in assemblies.keys():
                    assemblies[a] = {}
        else:
            assem_list = []
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
            for a in assem_list:
                if sample not in assemblies[a].keys():
                    assemblies[a][sample] = {unit: {}}
                if unit not in assemblies[a][sample].keys():
                    assemblies[a][sample][unit] = {}
                if r2:
                    assemblies[a][sample][unit]["R1"] = [R1_p]
                    assemblies[a][sample][unit]["R2"] = [R2_p]
                else:
                    assemblies[a][sample][unit]["se"] = [se_p]
        # If there is no 'fq2' column, add the single file path as 'se'
        else:
            samples[sample][unit]["se"] = R1
            for a in assem_list:
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
    if "uppmax.uu.se" in hostname:
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
    tools = [config["annotation"][key] for key in config["annotation"].keys()]
    assems = [config["assembly"]["metaspades"], config["assembly"]["megahit"]]
    config["run_assembly"] = False
    if True in tools and True in assems:
        config["run_annotation"] = True
        # if True also assume the user wants assembly
        config["run_assembly"] = True
        # Set megahit as default unless metaspades is set
        if not config["assembly"]["megahit"] and not config["assembly"]["metaspades"]:
            config["assembly"]["megahit"] = True
    return config


def check_assembly(config, assemblies):
    """
    Check assemblies and config settings

    :param config: Snakemake config
    :param assemblies: Assembly dictionary
    :return: Tuple of updated config and assembly dict
    """
    config["assembler"] = "Megahit" if config["megahit"] else "Metaspades"
    if len(assemblies) > 0:
        if config["assembly"]["metaspades"]:
            # Remove single-end only assemblies
            # that Metaspades won't be able to run
            assemblies = filter_metaspades_assemblies(assemblies)
    return assemblies


def check_classifiers(config):
    """
    Set paths and params specific to classifiers

    :param config: Snakemake config
    :return: Updated config dict
    """
    # Add read-based config info
    config["centrifuge"]["index_path"] = ""
    config["centrifuge"]["base"] = ""
    config["centrifuge"]["dir"] = ""
    if config["classification"]["centrifuge"]:
        # Check if custom database exists
        custom = expand("{b}.{i}.cf", b=config["centrifuge"]["custom"], i=[1, 2, 3])
        if list(set([os.path.exists(x) for x in custom]))[0]:
            config["centrifuge"]["index_path"] = config["centrifuge"]["custom"]
        # If not, use prebuilt default
        else:
            p = config["centrifuge"]["prebuilt"]
            config["centrifuge"]["index_path"] = "resources/centrifuge/" + p
        # Set centrifuge index config variables
        index_path = config["centrifuge"]["index_path"]
        config["centrifuge"]["dir"] = os.path.dirname(index_path)
        config["centrifuge"]["base"] = bn(index_path)

    config["kraken"]["index_path"] = ""
    config["kraken"]["mem"] = ""
    if config["classification"]["kraken"]:
        # Check if custom database exists
        custom = expand(
            config["kraken"]["custom"] + "/{n}.k2d", n=["hash", "opts", "taxo"]
        )
        if list(set(os.path.exists(x) for x in custom))[0]:
            config["kraken"]["index_path"] = config["kraken"]["custom"]
        # If not, use prebuilt or standard
        elif config["kraken"]["standard_db"]:
            config["kraken"]["index_path"] = "resources/kraken/standard"
        else:
            config["kraken"]["index_path"] = (
                "resources/kraken/prebuilt/" + config["kraken"]["prebuilt"]
            )
        if config["kraken"]["reduce_memory"]:
            config["kraken"]["mem"] += "--memory-mapping"
    return config


# preprocessing functions


def preprocessing_input(config):
    if config["run_preprocessing"] or config["preprocessing"]["fastqc"]:
        return config["paths"]["results"] + "/report/samples_report.html"
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


def get_all_files(samples, directory, suffix="", nested=False):
    """
    Returns a list of files based on samples and directory

    :param samples: Samples dictionary
    :param directory: Directory to find files in
    :param suffix: Suffix of files to return
    :param nested: If True look for files inside sample_run directory
    :return:
    """
    files = []
    if type(directory) == list:
        directory = directory[0]
    for sample in samples:
        for unit in samples[sample].keys():
            if nested:
                d = "{}/{}_{}".format(directory, sample, unit)
            else:
                d = "{}".format(directory)
            if is_pe(samples[sample][unit]):
                files.append(d + "/{}_{}_pe{}".format(sample, unit, suffix))
            else:
                files.append(d + "/{}_{}_se{}".format(sample, unit, suffix))
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
        files = expand(
            config["paths"]["results"]
            + "/intermediate/fastqc/{sample}_{unit}_{pair}{PREPROCESS}_fastqc.zip",
            sample=sample,
            unit=unit,
            pair=pairs,
            PREPROCESS=pre,
        )
        return files
    return []


def get_trim_logs(sample, unit, pairs, config, d):
    if (
        not config["preprocessing"]["trimmomatic"]
        and not config["preprocessing"]["cutadapt"]
    ):
        return []
    if config["preprocessing"]["trimmomatic"]:
        trimmer = "trimmomatic"
    else:
        trimmer = "cutadapt"
    files = expand(
        config["paths"]["results"]
        + "/intermediate/preprocess/{sample}_{unit}_{pair}{s}.{trimmer}.log",
        sample=sample,
        unit=unit,
        pair=pairs,
        s=d["trimming"],
        trimmer=trimmer,
    )
    return files


def get_filt_logs(sample, unit, seq_type, config, d):
    if not config["preprocessing"]["phix_filter"]:
        return []
    files = expand(
        config["paths"]["results"]
        + "/intermediate/preprocess/{sample}_{unit}_PHIX_{seq_type}{s}.log",
        sample=sample,
        unit=unit,
        seq_type=seq_type,
        s=d["phixfilt"],
    )
    return files


def get_sortmerna_logs(sample, unit, seq_type, config):
    if not config["preprocessing"]["sortmerna"]:
        return []
    files = expand(
        config["paths"]["results"]
        + "/intermediate/preprocess/{sample}_{unit}_{seq_type}.sortmerna.log",
        sample=sample,
        unit=unit,
        seq_type=seq_type,
    )
    if not config["sortmerna"]["remove_filtered"]:
        d = {"non_rRNA": "rRNA", "rRNA": "non_rRNA"}
        filtered = d[config["sortmerna"]["keep"]]
        if seq_type == "pe":
            files.append(
                f"{config['paths']['results']}/intermediate/preprocess/{sample}_{unit}.sortmerna_unmerge.{filtered}.log"
            )
        else:
            files.append(
                f"{config['paths']['results']}/intermediate/preprocess/{sample}_{unit}_se.sortmerna_zip_{filtered}.log"
            )
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
        trimsettings = " {} ILLUMINACLIP:{}:{}".format(
            pre_adapter_params, adapter, adapter_params
        )
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
    # Quit and warn if all assemblies have been removed
    if len(d) == 0:
        sys.exit(
            """
WARNING: Metaspades requires paired-end data but all specified assemblies
only contain single-end data. Exiting...
        """
        )
    return d


def get_all_assembly_files(assembly_dict):
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
            files.append(
                results_path
                + "/assembly/{}/mapping/{}_{}_{}{}.bam".format(
                    g, sample, unit, seq_type, POSTPROCESS
                )
            )
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

    for record in SeqIO.parse(gz.open(f, "rt"), "fastq"):
        record.id = "{}_{}".format(i, record.id)
        SeqIO.write(record, fh, "fastq")
    return fh


# binning functions


def binning_input(config, report=False):
    """
    Generates input list for the binning part of the workflow

    :param config: Snakemake config
    :param report: Whether to gather input for the bin report rule
    :return:
    """
    bin_input = []
    if len(get_binners(config)) > 0:
        bin_input.append(
            config["paths"]["results"] + "/report/binning/binning_summary.tsv"
        )
    if config["binning"]["checkm"]:
        bin_input.append(config["paths"]["results"] + "/report/checkm/checkm.stats.tsv")
        # Don't include profile in report
        if not report:
            bin_input.append(
                config["paths"]["results"] + "/report/checkm/checkm.profiles.tsv"
            )
    if config["binning"]["gtdbtk"]:
        bin_input.append(
            config["paths"]["results"] + "/report/gtdbtk/gtdbtk.summary.tsv"
        )
        bin_input.append(
            config["paths"]["results"] + "/report/bin_annotation/tRNA.total.tsv"
        )
        bin_input.append(
            config["paths"]["results"] + "/report/bin_annotation/rRNA.types.tsv"
        )
    config["fastani"]["ref_genomes"] = {}
    if config["binning"]["fastani"]:
        bin_input.append(
            config["paths"]["results"] + "/report/binning/genome_clusters.tsv"
        )
        # read list of genome references if path exists
        if os.path.exists(config["fastani"]["ref_list"]):
            _ = pd.read_csv(
                config["fastani"]["ref_list"],
                index_col=0,
                sep="\t",
                header=None,
                names=["genome_id", "url"],
            )
            # filter genome list
            _ = _.loc[(_["url"].str.contains("ftp")) | (_["url"].str.contains("http"))]
            config["fastani"]["ref_genomes"] = _.to_dict()["url"]
        else:
            config["fastani"]["ref_genomes"] = {}
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
                r = "R1"
            else:
                r = "se"
            f = config["paths"][
                "results"
            ] + "/intermediate/preprocess/{sample}_{unit}_{r}{p}.fastq.gz".format(
                sample=sample, unit=unit, r=r, p=p
            )
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
    if config["checkm"]["reduced_tree"]:
        return "-r"
    return ""


def get_binners(config):
    """
    Return a list of binners used
    :param config:
    :return:
    """
    binners = []
    if config["binning"]["metabat"]:
        binners.append("metabat")
    if config["binning"]["concoct"]:
        binners.append("concoct")
    if config["binning"]["maxbin"]:
        binners.append("maxbin")
    return binners


def get_fields(f, index):
    """
    Extract assembly, binner and length fields from file path
    :param f: Input file
    :param index: Path split index starting from right
    :return:
    """
    items = f.split("/")
    return items[index - 2], items[index - 1], items[index]


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


def concatenate(input, index):
    """
    Concatenate bin info dataframes
    :param input:
    :return:
    """
    df = pd.DataFrame()
    for i, f in enumerate(input):
        binner, group, l = get_fields(f, index)
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
    results = config["paths"]["results"]
    if not config["assembly"]["megahit"] and not config["assembly"]["metaspades"]:
        return input
    for assembly in assemblies.keys():
        # Add orfcalling results
        input.append(f"{results}/annotation/{assembly}/final_contigs.gff")
        if config["annotation"]["infernal"]:
            input.append(f"{results}/annotation/{assembly}/final_contigs.cmscan")
        if config["annotation"]["tRNAscan"]:
            input.append(f"{results}/annotation/{assembly}/tRNA.out")
        # Add EGGNOG annotation
        if config["annotation"]["eggnog"]:
            input += expand(
                "{results}/annotation/{assembly}/{db}.parsed.{norm_method}.tsv",
                results=[results],
                assembly=[assembly],
                db=["enzymes", "pathways", "kos", "modules"],
                norm_method=["counts", "rpkm", "TMM", "RLE", "CSS"],
            )
        # Add PFAM annotation
        if config["annotation"]["pfam"]:
            if config["annotation"]["splits"] > 1:
                input += expand(
                    "{results}/annotation/{assembly}/{assembly}.pfam.gathered",
                    results=[results],
                    assembly=[assembly],
                )
            input += expand(
                "{results}/annotation/{assembly}/pfam.parsed.{norm_method}.tsv",
                results=[results],
                assembly=[assembly],
                norm_method=["counts", "rpkm", "TMM", "RLE", "CSS"],
            )
        # Add taxonomic annotation
        if config["annotation"]["taxonomy"]:
            input += expand(
                "{results}/annotation/{assembly}/{seqTaxDB}.orfs.taxonomy.{counts_type}.tsv",
                results=[results],
                assembly=[assembly],
                seqTaxDB=config["taxonomy"]["database"],
                counts_type=["counts", "rpkm"],
            )
        # Add Resistance Gene Identifier output
        if config["annotation"]["rgi"]:
            input += expand(
                "{results}/annotation/{assembly}/rgi.parsed.{norm_method}.tsv",
                results=[results],
                assembly=[assembly],
                norm_method=["counts", "rpkm", "TMM", "RLE", "CSS"],
            )
    return input


def parse_cmout(f):
    with open(f) as fh:
        lines = []
        idnums = {}
        for i, line in enumerate(fh):
            if line[0] == "#":
                continue
            line = line.rstrip()
            line = line.rsplit()
            indices = [1, 2, 3, 5, 9, 10, 11, 14, 16, 17]
            (
                target_name,
                target_accession,
                query,
                clan,
                start,
                end,
                strand,
                gc,
                score,
                evalue,
            ) = [line[x] for x in indices]
            try:
                idnum = idnums[query]
            except KeyError:
                idnum = 0
                idnums[query] = idnum
            this_idnum = idnum + 1
            idnums[query] = this_idnum
            attributes = [
                "ID=" + query + "ncRNA_" + str(this_idnum),
                "Name=" + target_name,
                "Accession=" + target_accession,
                "Clan=" + clan,
                "GC=" + gc,
                "E-value=" + evalue,
            ]
            # seqid, source, type, start, end, score, strand, phase, attributes
            gffline = " ".join(
                [
                    query,
                    "cmscan",
                    "ncRNA",
                    start,
                    end,
                    score,
                    strand,
                    ".",
                    ";".join(attributes),
                ]
            )
            lines.append(gffline)
    return lines


# classification functions


def classify_input(config):
    f = []
    results = config["paths"]["results"]
    if config["classification"]["kraken"]:
        f.append(f"{results}/report/kraken/kraken.krona.html")
    if config["classification"]["metaphlan"]:
        f.append(f"{results}/report/metaphlan/metaphlan.html")
    if config["classification"]["centrifuge"]:
        f.append(f"{results}/report/centrifuge/centrifuge.krona.html")
    return f


def krona_input(config, samples, classifier):
    input_string = ""
    results = config["paths"]["results"]
    files = get_all_files(samples, f"{results}/{classifier}", ".kreport")
    for f in files:
        sample_unit = bn(f).replace("_pe.kreport", "").replace("_se.kreport", "")
        input_string += f" {f},{sample_unit}"
    return input_string


def metaphlan_krona_string(input):
    """
    Returns the input string (<file1>,<name1> ... <fileN>,<nameN>) for krona
    :param input: krona-ready metaphlan tables
    :return: krona-formatted input string
    """
    s = []
    for f in input:
        name = bn(f).replace("_pe.krona", "").replace("_se.krona", "")
        s.append(f"{f},{name}")
    return " ".join(s)
