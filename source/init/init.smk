import pandas as pd, subprocess
import platform
import glob, yaml, os
from os.path import join as opj
from source.utils.parse_samplelist import parse_samplelist, check_sequencing_type

# please note that the run identifier is supposed to be numbers only!
wildcard_constraints:
    run="\d+",
    pair="se|R[12]",
    seq_type="[sp]e"

def is_pe(d):
    if "R2" in d.keys():
        return True
    return False

def get_all_files(samples, dir, suffix=""):
    files = []
    for sample in samples:
        for run in samples[sample].keys():
            if is_pe(samples[sample][run]):
                files.append(opj(dir, "{}_{}_pe{}".format(sample,run,suffix)))
            else:
                files.append(opj(dir, "{}_{}_se{}".format(sample,run,suffix)))
    return files

###########################################################
# Store config parameters to write to pipeline_report.txt #
###########################################################
config_params = []

###############
# Check paths #
###############
pythonpath = sys.executable
envdir = '/'.join(pythonpath.split("/")[0:-2])
system = platform.system()
config["system"] = system
config_params.append((" - System",system))

config["tmpdir"] = config["temp_path"]

# Figure out trimmomatic home if configuration path not available
if not os.path.exists(config["trimmomatic_home"]) and config["trimmomatic"]:
    config["trimmomatic_home"] = envdir+"/share/trimmomatic/"
    if not os.path.exists(config["trimmomatic_home"]):
        print("Trimmomatic path not found",file=sys.stderr)

# Locate picard.jar if markduplicates is used
if not config["picard_jar"] and config["markduplicates"]:
    for line in shell("find {envdir} -name picard.jar", iterable = True):
        basename = os.path.basename(line)
        if basename == "picard.jar":
            config["picard_jar"] = line
            config["picard_path"] = os.path.dirname(line)

# Check that workflow is actually running on Uppmax if so specified
hostname = platform.node()
if 'uppmax.uu.se' in hostname:
    config["runOnUppMax"] = 'yes'
    config_params.append((" - Cluster node",hostname))
    # Set scratch_path to $TMPDIR and use os.path.expandvars(config["scratch_path"]) to call it from rules.
    config["scratch_path"] = "$TMPDIR"
    config["tmpdir"] = "$TMPDIR"
else:
    config["runOnUppMax"] = "no"
config_params.append((" - Temporary path",os.path.abspath(config["tmpdir"])))
config_params.append((" - Scratch path",os.path.abspath(config["scratch_path"])))
config_params.append((" - Intermediate path", os.path.abspath(config["intermediate_path"])))
config_params.append((" - Results path",os.path.abspath(config["results_path"])))
config_params.append((" - Resource path",os.path.abspath(config["resource_path"])))

# Check whether to set annotation downstream of assembly
if config["pfam"] or config["taxonomic_annotation"] or config["infernal"] or config["eggnog"]:
    config["annotation"] = True
    #if True also assume the user wants assembly
    config["assembly"] = True
else:
    config["annotation"] = False

# Check that taxonomic dbtype is correctly entered
if config["taxonomic_annotation"]:
    if config["diamond_dbtype"] not in ["nr","uniref50","uniref90","uniref100"]:
        config["taxonomic_annotation"] = False
        config["diamond_protdb"] = ""
    else:
        config["diamond_protdb"] = opj(config["diamond_dbpath"],"diamond_"+config["diamond_dbtype"]+".dmnd")
        taxonmap = opj(config["taxdb"],"{}.accession2taxid.gz".format(config["diamond_dbtype"]))
else:
    config["diamond_protdb"] = ""

#######################################################
# Figure out pre- and post-processing to be performed #
#######################################################
PREPROCESS=""
POSTPROCESS=""
preprocess_suffices = {"sortmerna": "", "trimming": "", "phixfilt": "", "fastuniq": ""}

# SortMeRNA rRNA filtering
if config["sortmerna"]:
    PREPROCESS+=".sortmerna"
    if config["sortmerna_keep"].lower() in ["non_rrna", "rrna"]:
        preprocess_suffices["trimming"] = ".sortmerna"
    config_params.append((" - Preprocessing","SortMeRNA"))
    config_params.append(("   - keep", "{}".format(config["sortmerna_keep"])))
    config_params.append(("   - rRNA databases","{}".format(" ".join(config["sortmerna_dbs"]))))

# Trimming
if config["trimmomatic"]:
    PREPROCESS+=".trimmomatic"
    preprocess_suffices["phixfilt"] = preprocess_suffices["trimming"]+".trimmomatic"
    config_params.append((" - Preprocessing","Trimmomatic"))
    for key in ["pe_adapter_params", "se_adapter_params", "pe_pre_adapter_params", "pe_post_adapter_params", "se_pre_adapter_params", "se_post_adapter_params"]:
        config_params.append(("    - ","{}: {}".format(key,config[key])))
elif config["cutadapt"]:
    PREPROCESS+=".cutadapt"
    preprocess_suffices["phixfilt"] = preprocess_suffices["trimming"]+".cutadapt"
    config_params.append((" - Preprocessing","Cutadapt"))
    for key in ["adapter_sequence", "rev_adapter_sequence"]:
        config_params.append(("    - ","{}: {}".format(key,config["cutadapt"])))
else:
    preprocess_suffices["phixfilt"] = preprocess_suffices["trimming"]

# Filtering
if config["phix_filter"]:
    preprocess_suffices["fastuniq"] = preprocess_suffices["phixfilt"]+".phixfilt"
    PREPROCESS+=".phixfilt"
    config_params.append((" - Preprocessing","PhiX filtering"))
else:
    preprocess_suffices["fastuniq"] = preprocess_suffices["phixfilt"]

# Deduplication
if config["fastuniq"]:
    PREPROCESS+=".fastuniq"
    config_params.append((" - Preprocessing","Fastquniq"))

if PREPROCESS!="":
    config["preprocess"] = True
else:
    config["preprocess"] = False

# Duplicate removal using picard and MarkDuplicates
if config["markduplicates"]:
    POSTPROCESS += ".markdup"
    config_params.append((" - Alignment postprocessing","Markduplicates"))

#####################
# Parse sample list #
#####################
df = None
if(os.path.isfile(config["sample_list"])):
    samples,assemblyGroups = parse_samplelist(config["sample_list"],config,PREPROCESS)
    config_params.append((" - Sample list",os.path.abspath(config["sample_list"])))
    # Check that sequencing_type matches with all files in the input if
    seq_type = check_sequencing_type(samples)
    config["seq_type"] = seq_type
    config_params.append((" - Data Type (paired, single, mixed)",seq_type))

    if len(assemblyGroups) > 0 and config["assembly"]:
        config_params.append((" - Assemblies to generate", len(assemblyGroups)))
        config_params.append(("   - Assembler", "megahit"))
        config_params.append(("   - Keep intermediate contigs", config["megahit_keep_intermediate"]))
        config_params.append(("   - Assembly additional params", config["megahit_additional_settings"]))
        # Add information on binning
        if config["maxbin"]:
            config_params.append((" - Genome binning", "MaxBin2"))
        if config["concoct"]:
            config_params.append((" - Genome binning", "CONCOCT"))
            config_params.append(("   - Min contig length",",".join(str(x) for x in config["min_contig_length"])))
else:
    print("Could not read the sample list file, wont be able to run the pipeline, tried "+config["sample_list"])
    samples = {}
    assemblyGroups = {}
    mapping = {}

# Add config information related to orf-calling
if config["infernal"]:
    config_params.append((" - Non coding RNA database", os.path.abspath(opj(config["infernal_dbpath"],"Rfam.cm"))))

# Add config information related to annotation
if config["pfam"]:
    config_params.append((" - PFAM database directory", os.path.abspath(opj(config["resource_path"],"pfam"))))
if config["taxonomic_annotation"]:
    config_params.append((" - Diamond database", os.path.abspath(opj(config["diamond_dbpath"], "{}.fasta".format(config["diamond_dbtype"])))))

# Add read-based config info
if config["centrifuge"]:
    # Check if custom database exists
    custom = expand("{b}.{i}.cf", b=config["centrifuge_custom"], i=[1,2,3])
    if list(set([os.path.exists(x) for x in custom]))[0]:
        config["centrifuge_index_path"] = config["centrifuge_custom"]
    # If not, use prebuilt default
    else:
        config["centrifuge_index_path"] = "resources/classify_db/centrifuge/{}".format(config["centrifuge_prebuilt"])
    config_params.append((" - Read classifier","Centrifuge"))
    # Set centrifuge index config variables
    config['centrifuge_dir'] = os.path.dirname(config['centrifuge_index_path'])
    config['centrifuge_base'] = os.path.basename(config['centrifuge_index_path'])
if config["kaiju"]:
    config_params.append((" - Read classifier","Kaiju"))
if config["kraken"]:
    # Check if custom database exists
    custom = expand(opj(config["kraken_custom"], "{n}.k2d"), n=["hash","opts","taxo"])
    if list(set(os.path.exists(x) for x in custom))[0]:
        config["kraken_index_path"] = config["kraken_custom"]
    # If not, use prebuilt default
    else:
        config["kraken_index_path"] = "resources/classify_db/{}".format(config["kraken_prebuilt"])
    if config["kraken_reduce_memory"]:
        config["kraken_params"] = "--memory-mapping"
    else:
        config["kraken_params"] = ""
    config_params.append((" - Read classifier","Kraken"))
if config["reference_map"]:
    config_params.append((" - Reference based mapping","True"))


config_params.append((" - Configfiles", "{}".format(",".join(workflow.configfiles))))
localrules: write_config
# Rule for writing pipeline configuration to file
rule write_config:
    output:
        config["pipeline_config_file"]
    run:
        import yaml
        dirname = config["results_path"]
        shell("mkdir -p {dirname}")
        with open(output[0], 'w') as fh_out:
            fh_out.write("{}\n{}\n".format("Pipeline information","-"*len("Pipeline information")))
            for key,val in config_params:
                fh_out.write("{}: {}\n".format(key,val))
            fh_out.write(("\n{}\n{}\n".format("Configuration file settings","-"*len("Configuration file settings"))))
            yaml.safe_dump(config, fh_out)
            fh_out.write(("\n{}\n{}\n".format("Software versions","-"*len("Software versions"))))
            for line in shell("conda list", iterable = True):
                if line[0] == "#":
                    continue
                items = line.rstrip().rsplit()
                if len(items) == 3:
                    items.append("")
                fh_out.write(" - {}: version {}, build {} {}\n".format(items[0], items[1], items[2], items[3]))


rule download_examples:
    output:
        expand("examples/data/{sample}_100000_R{i}.fastq.gz",
            sample = ["anterior_nares", "buccal_mucosa", "retr_crease", "stool"], i = ["1","2"])
    params:
        url = "https://bitbucket.org/johnne/metagenomic-mocks/raw/995102fbc6e19963f9b4861da0050b41c04a0063/results"
    run:
        for f in output:
            fn=os.path.basename(f)
            shell("curl -L -o examples/data/{fn} {params.url}/{fn}")