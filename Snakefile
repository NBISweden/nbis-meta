from snakemake.utils import min_version, validate
from snakemake.exceptions import WorkflowError
min_version("4.4.0")

# Snakemake workflow for various types of metagenomics analyses.
# See documentation at https://bitbucket.org/scilifelab-lts/nbis-meta

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
    return

shell.prefix("")
configfile: "config.yaml"
try:
    validate(config, "config/config.schema.yaml")
except WorkflowError as e:
    parse_validation_error(e)
    sys.exit()
workdir: config["workdir"]

# First load init file to set up samples and variables
include: "source/init/init.smk"
pipeline_report = config["pipeline_config_file"]
###########
# Targets #
###########
inputs = [pipeline_report]
db_done = opj(config["results_path"],"progress","db.done")
preprocess_done = opj(config["results_path"],"progress","preprocess.done")
assembly_done = opj(config["results_path"],"progress","assembly.done")
binning_done = opj(config["results_path"],"progress","binning.done")
checkm_db_done = opj(config["results_path"],"progress","checkm_db.done")
annotation_done = opj(config["results_path"],"progress","annotation.done")
kraken_db_done = opj(config["results_path"],"progress","kraken_db.done")
kraken_classify_done = opj(config["results_path"],"progress","kraken_classify.done")
metaphlan2_db_done = opj(config["results_path"],"progress","metaphlan2_db.done")
metaphlan2_classify_done = opj(config["results_path"],"progress","metaphlan2_classify.done")
kaiju_db_done = opj(config["results_path"],"progress","kaiju_db.done")
kaiju_classify_done = opj(config["results_path"],"progress","kaiju_classify.done")
centrifuge_db_done = opj(config["results_path"],"progress","centrifuge_db.done")
refmap_done = opj(config["results_path"],"progress","refmap.done")
centrifuge_classify_done = opj(config["results_path"],"progress","centrifuge_classify.done")

# Download and format databases for annotation
include: "source/workflow/DB"
inputs.append(db_done)

# Preprocess raw input (if no preprocessing, just produce the sample report for raw data)
include: "source/workflow/Preprocess"
inputs.append(preprocess_done)

# Assemble
if config["assembly"]:
    include: "source/workflow/Assembly"
    inputs.append(assembly_done)
    # Rule sets that depend on de-novo assembly
    # Annotate
    if config["annotation"]:
        include: "source/workflow/Annotation"
        inputs.append(annotation_done)
    # Binning
    if config["binning"]:
        inputs.append(binning_done)
        include: "source/workflow/Binning"
      # Checkm
    if config["checkm"]:
        inputs.append(checkm_db_done)
        include: "source/workflow/CheckmDB"
    else:
        config["checkm"] = False
# Kraken
if config["kraken"]:
    # Download and process kraken datatbase
    include: "source/workflow/KrakenDB"
    # Kraken classify samples
    include: "source/workflow/KrakenClassify"
    inputs += [kraken_db_done, kraken_classify_done]
# Kaiju
if config["kaiju"]:
    # Process the Kaiju database
    include: "source/workflow/KaijuDB"
    # Kaiju classify samples
    include: "source/workflow/KaijuClassify"
    inputs += [kaiju_db_done, kaiju_classify_done]
# Metaphlan2
if config["metaphlan2"]:
    include: "source/workflow/Metaphlan2DB"
    include: "source/workflow/Metaphlan2Classify"
    inputs += [metaphlan2_db_done, metaphlan2_classify_done]
# Centrifuge
if config["centrifuge"]:
    include: "source/workflow/CentrifugeDB"
    include: "source/workflow/CentrifugeClassify"
    inputs += [centrifuge_db_done, centrifuge_classify_done]

# Reference-based mapping
if config["reference_map"]:
    # Use centrifuge to download genomes for reference mapping
    # So set run_centrifuge to True
    config["centrifuge"] = True
    include: "source/workflow/Map"
    inputs.append(refmap_done)

# master target rule
rule all:
    input: inputs

# db target rule
rule db:
  input: db_done

# preprocess target rule
rule preprocess:
    input: pipeline_report, preprocess_done

# assembly target rule
rule assembly:
    input: pipeline_report, preprocess_done, assembly_done

# annotation target rule
rule annotation:
    input: pipeline_report, preprocess_done, db_done, assembly_done, annotation_done

# centrifuge
rule centrifuge_db:
    input: centrifuge_db_done
rule centrifuge_classify:
    input: pipeline_report, preprocess_done, centrifuge_classify_done

# kraken
rule kraken_db:
    input: kraken_db_done
rule kraken_classify:
    input: pipeline_report, preprocess_done, kraken_classify_done

# kaiju
rule kaiju_db:
    input: kaiju_db_done
rule kaiju_classify:
    input: pipeline_report, preprocess_done, kaiju_db_done, kaiju_classify_done

# metaphlan2
rule metaphlan2_db:
    input: metaphlan2_db_done
rule metaphlan2_classify:
    input: pipeline_report, preprocess_done, metaphlan2_db_done, metaphlan2_classify_done

# binning
rule binning:
    input: pipeline_report, preprocess_done, assembly_done, binning_done

# checkm
rule checkm_db:
    input: pipeline_report, checkm_db_done

# Reference based database
rule refmap:
    input: pipeline_report, centrifuge_db_done, preprocess_done, refmap_done
