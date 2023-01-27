from snakemake.utils import validate
import pandas as pd
import platform
import os


include: "../scripts/common.py"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3:4.8.2"


##### load and validate config #####

if os.path.exists("config/config.yaml"):

    configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml", set_default=True)


# Set results and temp directories
results = config["paths"]["results"]
temppath = config["paths"]["temp"]

##### generate preprocess/postprocess strings #####
PREPROCESS, POSTPROCESS, preprocess_suffices, config = prepost_string(config)

##### load and validate samples #####
df = pd.read_csv(config["sample_list"], sep="\t")
validate(df, schema="../schemas/samples.schema.yaml")

##### parse samples #####
samples, assemblies = parse_samples(df, config, PREPROCESS)

##### workflow settings #####


wildcard_constraints:
    unit="\d+",
    pair="se|R[12]",
    seq_type="[sp]e",
    binner="[a-z]+",
    group="\w+",
    l="\d+",
    counts_type="(counts|rpkm)",
    norm_method="(TMM|RLE)",


config = check_uppmax(config)
config = check_annotation(config)
assemblies = check_assembly(config, assemblies)
config = check_classifiers(config)
