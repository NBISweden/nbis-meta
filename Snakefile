# Snakemake workflow for various types of metagenomics analyses.
include: "source/rules/common.smk"
include: "source/rules/examples.smk"

#############
## Targets ##
#############
inputs = []
# Download and format databases for annotation
db_input = []
include: "source/workflow/db.smk"
inputs += db_input
# Preprocess raw input (if no preprocessing, just produce the sample report for raw data)
preprocess_input = []
include: "source/workflow/preprocess.smk"
inputs += preprocess_input
# Assemble
assembly_input = []
annotation_input = []
binning_input = []
if config["assembly"]:
    include: "source/workflow/assembly.smk"
    inputs += assembly_input
    # Rule sets that depend on de-novo assembly
    # Annotate
    if config["annotation"]:
        include: "source/workflow/annotation.smk"
        inputs += annotation_input
    # Binning
    if config["maxbin"] or config["concoct"] or config["metabat"]:
        include: "source/workflow/binning.smk"
        inputs += binning_input
# Kraken
kraken_input = []
kraken_db_input = []
if config["kraken"]:
    # Download and process kraken datatbase
    include: "source/workflow/kraken_db.smk"
    # Kraken classify samples
    include: "source/workflow/kraken.smk"
    inputs += kraken_input + kraken_db_input
# Metaphlan2
metaphlan_input = []
metaphlan_db_input = []
if config["metaphlan"]:
    include: "source/workflow/metaphlan_db.smk"
    include: "source/workflow/metaphlan.smk"
    inputs += metaphlan_input + metaphlan_db_input
# Centrifuge
centrifuge_input = []
centrifuge_db_input = []
if config["centrifuge"]:
    include: "source/workflow/centrifuge_db.smk"
    include: "source/workflow/centrifuge.smk"
    inputs += centrifuge_input + centrifuge_db_input
if config["centrifuge"] or config["kraken"] or config["metaphlan"]:
    include: "source/rules/Classify/krona.smk"

###########
## RULES ##
###########
# master target rule
rule all:
    input: inputs

# db target rule
rule db:
  input: db_input

# preprocess target rule
rule preprocess:
    input: preprocess_input

# assembly target rule
rule assembly:
    input: preprocess_input, assembly_input

# annotation target rule
rule annotation:
    input: preprocess_input, db_input, assembly_input, annotation_input

# centrifuge
rule centrifuge_db:
    input: centrifuge_db_input
rule centrifuge_classify:
    input: preprocess_input, centrifuge_input

# kraken
rule kraken_db:
    input: kraken_db_input
rule kraken_classify:
    input: preprocess_input, kraken_input

# metaphlan2
rule metaphlan_db:
    input: metaphlan_db_input
rule metaphlan_classify:
    input: preprocess_input, metaphlan_db_input, metaphlan_input

# binning
rule binning:
    input: preprocess_input, assembly_input, binning_input