from snakemake.utils import validate
from snakemake.exceptions import WorkflowError
import gzip as gz
from Bio import SeqIO
import pandas as pd
import platform
import os
from os.path import join as opj
from source.utils.parse_samplelist import parse_samplelist, check_sequencing_type


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


wildcard_constraints:
    run="\d+",
    pair="se|R[12]",
    seq_type="[sp]e"

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
    print("Could not read the sample list file, wont be able to run the pipeline, tried "+config["sample_list"])
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
    config_params.append((" - Read classifier","Centrifuge"))
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

rule download_synthetic:
    """
    Download pre-made synthetic metagenome from Zenodo
    """
    output:
        R1 = temp("examples/data/synthetic_1.fastq.gz"),
        R2 = temp("examples/data/synthetic_2.fastq.gz")
    params:
        tar = "examples/data/synthetic.tar.gz",
        url = "https://zenodo.org/record/3737112/files/synthetic.tar.gz?download=1"
    conda:
        "../../envs/examples.yml"
    shell:
         """
         curl -L -s -o {params.tar} {params.url}
         tar -C examples/data/ -xf {params.tar}
         rm {params.tar}
         """

rule generate_examples:
    """
    Use seqtk to subsample the synthetic metagenome into examples
    """
    input:
        "examples/data/synthetic_{i}.fastq.gz"
    output:
        "examples/data/{sample}_{s}_R{i}.fastq.gz"
    conda:
        "../../envs/examples.yml"
    shell:
         """
         seqtk sample -s {wildcards.s} {input[0]} 100000 | gzip -c > {output[0]}
         """