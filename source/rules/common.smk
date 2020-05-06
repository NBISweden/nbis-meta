from snakemake.utils import validate
from snakemake.exceptions import WorkflowError
import gzip as gz
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import platform
import os
from os.path import join as opj
from source.utils.parse_samplelist import parse_samplelist, check_sequencing_type


def link(target,link_name):
    target_abs = os.path.abspath(target)
    link_abs = os.path.abspath(link_name)
    shell("ln -s {target_abs} {link_abs}")


def get_interleaved(sample,runID):
    files = []
    if "interleaved" in samples[sample][runID].keys():
        inter = samples[sample][runID]["interleaved"]
        R1 = samples[sample][runID]["R1"]
        R2 = samples[sample][runID]["R2"]
        files.append(inter)
    else:
        files.append("")
    return files


def is_pe(d):
    if "R2" in d.keys():
        return True
    return False


def parse_validation_error(e):
    instance = ""
    print("ERROR VALIDATING CONFIG FILE")
    for item in str(e).split("\n"):
        item = item.replace('"','')
        if "ValidationError:" in item:
            print(item)
        if item[0:11] == "On instance":
            instance = item.replace("On instance['", "INCORRECT CONFIG AT: ")
            instance = instance.replace("']:","")
            print(instance)
    sys.exit()


def get_all_files(samples, dir, suffix="", nested=False):
    """
    Returns a list of files based on samples and directory

    :param samples: Samples dictionary
    :param dir: Directory to find files in
    :param suffix: Suffix of files to return
    :param nested: If True look for files inside sample_run directory
    :return:
    """
    files=[]
    for sample in samples:
        for run in samples[sample].keys():
            if nested:
                d = "{}/{}_{}".format(dir,sample,run)
            else:
                d = "{}".format(dir)
            if is_pe(samples[sample][run]):
                files.append(opj(d, "{}_{}_pe{}".format(sample,run,suffix)))
            else:
                files.append(opj(d, "{}_{}_se{}".format(sample,run,suffix)))
    return files


def get_fastqc_files(wildcards):
    """Get all fastqc output"""
    files = []
    for sample in samples.keys():
        for run in samples[sample].keys():
            for pair in samples[sample][run].keys():
                if pair in ["R1","R2","se"]:
                    if config["preprocess"]:
                        files.append(opj(config["intermediate_path"],
                            "fastqc","{}_{}_{}{}_fastqc.zip".format(sample,
                            run,pair,PREPROCESS)))
                    else:
                        files.append(opj(config["intermediate_path"],
                            "fastqc","{}_{}_{}_fastqc.zip".format(sample,
                            run,pair)))
    return files

def get_trim_logs(wildcards):
    """
    Get all trimming logs from Trimmomatic and/or cutadapt

    :param wildcards: wildcards from snakemake
    :return: list of files
    """
    files = []
    if not config["trimmomatic"] and not config["cutadapt"]:
        return files
    if config["trimmomatic"]:
        trimmer = "trimmomatic"
    elif config["cutadapt"]:
        trimmer = "cutadapt"
    for sample in samples.keys():
        for run in samples[sample].keys():
            for pair in samples[sample][run].keys():
                if pair in ["R1","R2","se"]:
                    logfile = opj(config["intermediate_path"],"preprocess",
                                  "{}_{}_{}{}.{}.log".format(sample,
                        run,pair,preprocess_suffices["trimming"],trimmer))
                    files.append(logfile)
    return files

def get_filt_logs(wildcards):
    """
    Get all filter logs from Phix filtering

    :param wildcards: wildcards from snakemake
    :return: list of files
    """
    files = []
    if not config["phix_filter"]: return files
    for sample in samples.keys():
        for run in samples[sample].keys():
            if "R2" in samples[sample][run].keys():
                logfile = opj(config["intermediate_path"],"preprocess",
                    "{}_{}_PHIX_pe{}.log".format(sample,run,
                                            preprocess_suffices["phixfilt"]))
            else:
                logfile = opj(config["intermediate_path"],"preprocess",
                    "{}_{}_PHIX_se{}.log".format(sample,run,
                                            preprocess_suffices["phixfilt"]))
            files.append(logfile)
    return files

def get_sortmerna_logs(wildcards):
    """
    Get all logs from SortMeRNA

    :param wildcards: wildcards from snakemake
    :return: list of files
    """
    files = []
    if not config["sortmerna"]:
        return files
    for sample in samples.keys():
        for run in samples[sample].keys():
            if "R2" in samples[sample][run].keys():
                logfile = opj(config["intermediate_path"],"preprocess",
                              "{}_{}_pe.sortmerna.log".format(sample,run))
            else:
                logfile = opj(config["intermediate_path"],"preprocess",
                              "{}_{}_se.sortmerna.log".format(sample,run))
            files.append(logfile)
    return files


def get_trimmomatic_string(seq_type):
    """
    Generates trimsetting string for Trimmomatic

    :param seq_type: PE or SE depending on sequencing type
    :return: trimsettings string
    """
    trim_adapters=config["trim_adapters"]
    adapter_fasta_dir="$CONDA_PREFIX/share/trimmomatic/adapters"
    adapter="{}/{}.fa".format(adapter_fasta_dir,
                            config["trimmomatic_{}_adapter".format(seq_type)])
    adapter_params=config["{}_adapter_params".format(seq_type)]
    pre_adapter_params=config["{}_pre_adapter_params".format(seq_type)]
    post_adapter_params=config["{}_post_adapter_params".format(seq_type)]
    trimsettings=pre_adapter_params
    if trim_adapters:
        trimsettings+=" ILLUMINACLIP:"+adapter+":"+adapter_params
    trimsettings+=" "+post_adapter_params
    return trimsettings


def get_sortmerna_ref_string(path, s):
    """
    Constructs the SortMeRNA --ref string

    :param path: Resource path from config
    :param s: Sortmerna databases from config
    :return: STRING,STRING formatted string
    """
    files=["{p}/rRNA_databases/{db}".format(p=path, db=db)
           for db in config["sortmerna_dbs"]]
    ref_string =":".join(
        ["{},{}".format(f,f) for f in files])
    return ref_string

## Assembly functions

def get_all_group_files(g):
  files=[]
  for sample in assemblyGroups[g].keys():
    for run in assemblyGroups[g][sample].keys():
      for pair in assemblyGroups[g][sample][run].keys():
        files.append(assemblyGroups[g][sample][run][pair][0])
  return files

def get_bamfiles(g):
  files=[]
  for sample in assemblyGroups[g].keys():
    for run in assemblyGroups[g][sample].keys():
      if "R2" in assemblyGroups[g][sample][run].keys():
        files.append(opj(config["results_path"],"assembly",g,"mapping",sample+"_"+run+"_pe"+POSTPROCESS+".bam"))
      else:
        files.append(opj(config["results_path"],"assembly",g,"mapping",sample+"_"+run+"_se"+POSTPROCESS+".bam"))
  return files

def rename_records(f, fh, i):
    """
    Prepends a number to read ids

    :param f: Input fastq file (gzipped)
    :param fh: Output filehandle
    :param i: File index to prepend to reads
    :return: Output filehandle
    """
    for record in SeqIO.parse(gz.open(f, 'rt'), 'fastq'):
        record.id="{}_{}".format(i, record.id)
        SeqIO.write(record, fh, "fastq")
    return fh

## Annotation functions

def parse_cmout(f):
    with open(f) as fh:
        lines = []
        idnums = {}
        for i, line in enumerate(fh):
            if line[0] == "#": continue
            line = line.rstrip()
            line = line.rsplit()
            indices = [1,2,3,5,9,10,11,14,16,17]
            target_name, target_accession, query, clan, start, end, strand, gc, score, evalue = [line[x] for x in indices]
            try:
                idnum = idnums[query]
            except KeyError:
                idnum = 0
                idnums[query] = idnum
            this_idnum = idnum+1
            idnums[query] = this_idnum
            attributes = ["ID="+query+"ncRNA_"+str(this_idnum),"Name="+target_name,"Accession="+target_accession,"Clan="+clan,"GC="+gc,"E-value="+evalue]
            # seqid, source, type, start, end, score, strand, phase, attributes
            gffline = " ".join([query,"cmscan","ncRNA",start,end,score,strand,".",";".join(attributes)])
            lines.append(gffline)
    return lines


## Quantification functions

def get_fc_files(wildcards, file_type):
    g=wildcards.group
    files=[]
    for sample in assemblyGroups[g].keys():
        for run in assemblyGroups[g][sample].keys():
            if "se" in assemblyGroups[g][sample][run].keys():
                files.append(opj(config["results_path"],"assembly",g,"mapping",sample+"_"+run+"_se.fc.{}.tab".format(file_type)))
            else:
                files.append(opj(config["results_path"],"assembly",g,"mapping",sample+"_"+run+"_pe.fc.{}.tab".format(file_type)))
    return files


def concat_files(files, gff_df):
    df=pd.DataFrame()
    for f in files:
        _df=pd.read_csv(f, index_col=0, sep="\t")
        df=pd.concat([df,_df], axis=1)
    df=pd.merge(df, gff_df, left_index=True, right_on="gene_id")
    df.drop("gene_id", axis=1, inplace=True)
    df.set_index("orf", inplace=True)
    return df


## WORKFLOW SETUP ##

wildcard_constraints:
    run="\d+",
    pair="se|R[12]",
    seq_type="[sp]e",
    group="\w+",
    l="\d+"

# Validate and get default settings from schema
try:
    validate(config, "../../config/config.schema.yaml")
except WorkflowError as e:
    parse_validation_error(e)

workdir: config["workdir"]

config["tmpdir"]=config["temp_path"]

# Check that workflow is actually running on Uppmax if so specified
hostname=platform.node()
if 'uppmax.uu.se' in hostname:
    config["runOnUppMax"]='yes'
    # Set scratch_path to $TMPDIR and use os.path.expandvars(config["scratch_path"]) to call it from rules.
    config["scratch_path"]="$TMPDIR"
    config["tmpdir"]="$TMPDIR"
else:
    config["runOnUppMax"]="no"

# Check whether to set annotation downstream of assembly
if config["pfam"] or config["taxonomic_annotation"] or config["infernal"] or config["eggnog"] or config["rgi"]:
    config["annotation"]=True
    #if True also assume the user wants assembly
    config["assembly"]=True
else:
    config["annotation"]=False

#######################################################
# Figure out pre- and post-processing to be performed #
#######################################################
PREPROCESS=""
POSTPROCESS=""
preprocess_suffices={"sortmerna": "", "trimming": "", "phixfilt": "", "fastuniq": ""}

# SortMeRNA rRNA filtering
if config["sortmerna"]:
    PREPROCESS+=".sortmerna"
    if config["sortmerna_keep"].lower() in ["non_rrna", "rrna"]:
        preprocess_suffices["trimming"]=".sortmerna"

# Trimming
if config["trimmomatic"]:
    PREPROCESS+=".trimmomatic"
    preprocess_suffices["phixfilt"]=preprocess_suffices["trimming"]+".trimmomatic"
elif config["cutadapt"]:
    PREPROCESS+=".cutadapt"
    preprocess_suffices["phixfilt"]=preprocess_suffices["trimming"]+".cutadapt"
else:
    preprocess_suffices["phixfilt"]=preprocess_suffices["trimming"]

# Filtering
if config["phix_filter"]:
    preprocess_suffices["fastuniq"]=preprocess_suffices["phixfilt"]+".phixfilt"
    PREPROCESS+=".phixfilt"
else:
    preprocess_suffices["fastuniq"]=preprocess_suffices["phixfilt"]

# Deduplication
if config["fastuniq"]:
    PREPROCESS+=".fastuniq"

if PREPROCESS!="":
    config["preprocess"]=True
else:
    config["preprocess"]=False

# Duplicate removal using picard and MarkDuplicates
if config["markduplicates"]:
    POSTPROCESS += ".markdup"


if os.path.isfile(config["sample_list"]):
    samples, assemblyGroups = parse_samplelist(config["sample_list"], config,
                                               PREPROCESS)
    # Check that sequencing_type matches with all files in the input if
    seq_type = check_sequencing_type(samples)
    config["seq_type"] = seq_type

    if len(assemblyGroups) > 0 and config["assembly"]:
        if config["metaspades"]:
            assembler = "metaspades"
        else:
            assembler = "megahit"
        # Add information on binning
        binning=False
        if config["maxbin"] or config["concoct"] or config["metabat"]:
            binning=True
else:
    print("Could not read the sample list file, wont be able to run the "
          "pipeline, tried {}".format(config["sample_list"]))
    samples={}
    assemblyGroups={}
    mapping={}

# Add read-based config info
if config["centrifuge"]:
    # Check if custom database exists
    custom=expand("{b}.{i}.cf", b=config["centrifuge_custom"], i=[1,2,3])
    if list(set([os.path.exists(x) for x in custom]))[0]:
        config["centrifuge_index_path"]=config["centrifuge_custom"]
    # If not, use prebuilt default
    else:
        config["centrifuge_index_path"]="resources/classify_db/centrifuge/{}".format(config["centrifuge_prebuilt"])
    # Set centrifuge index config variables
    config['centrifuge_dir']=os.path.dirname(config['centrifuge_index_path'])
    config['centrifuge_base']=os.path.basename(config['centrifuge_index_path'])

if config["kraken"]:
    # Check if custom database exists
    custom=expand(opj(config["kraken_custom"], "{n}.k2d"), n=["hash","opts","taxo"])
    if list(set(os.path.exists(x) for x in custom))[0]:
        config["kraken_index_path"]=config["kraken_custom"]
    # If not, use prebuilt default
    else:
        config["kraken_index_path"]="resources/classify_db/{}".format(config["kraken_prebuilt"])
    if config["kraken_reduce_memory"]:
        config["kraken_params"]="--memory-mapping"
    else:
        config["kraken_params"]=""