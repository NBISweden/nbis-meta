from snakemake.utils import validate
import pandas as pd
import platform
import os
from os.path import join as opj

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3:4.8.2"

##### load and validate config #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##### generate preprocess/postprocess strings #####
from scripts.common import prepost_string
PREPROCESS, POSTPROCESS, preprocess_suffices = prepost_string(config)

##### load and validate samples #####
df = pd.read_csv(config["sample_list"], sep="\t")
validate(df, schema="../schemas/samples.schema.yaml")

##### parse samples #####
from scripts.common import parse_samples
samples, assemblies = parse_samples(df, config, PREPROCESS)

##### workflow settings #####

wildcard_constraints:
    run="\d+",
    pair="se|R[12]",
    seq_type="[sp]e",
    group="\w+",
    l="\d+"

config["tmpdir"] = config["temp_path"]

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
    config["annotation"] = True
    #if True also assume the user wants assembly
    config["assembly"] = True
    if not config["megahit"] and not config["metaspades"]:
        config["megahit"] = True
else:
    config["annotation"]=False

# Check that sequencing_type matches with all files in the input if
if len(assemblies) > 0:
    if config["metaspades"]:
        from scripts.common import filter_metaspades_assemblies
        # Check that there are no single-end only assemblies
        # that Metaspades won't be able to run
        assemblies = filter_metaspades_assemblies(assemblies)
    # Add information on binning
    binning=False
    if config["maxbin"] or config["concoct"] or config["metabat"]:
        binning=True
    if not binning:
        config["checkm"] = False
        config["gtdbtk"] = False

# Add read-based config info
if config["centrifuge"]:
    # Check if custom database exists
    custom=expand("{b}.{i}.cf", b=config["centrifuge_custom"], i=[1,2,3])
    if list(set([os.path.exists(x) for x in custom]))[0]:
        config["centrifuge_index_path"]=config["centrifuge_custom"]
    # If not, use prebuilt default
    else:
        config["centrifuge_index_path"]="resources/centrifuge/{}".format(config["centrifuge_prebuilt"])
    # Set centrifuge index config variables
    config['centrifuge_dir']=os.path.dirname(config['centrifuge_index_path'])
    config['centrifuge_base']=os.path.basename(config['centrifuge_index_path'])

if config["kraken"]:
    # Check if custom database exists
    custom=expand(opj(config["kraken_custom"], "{n}.k2d"), n=["hash","opts","taxo"])
    if list(set(os.path.exists(x) for x in custom))[0]:
        config["kraken_index_path"]=config["kraken_custom"]
    # If not, use prebuilt or standard
    elif config["kraken_standard_db"]:
        config["kraken_index_path"] = opj(config["resource_path"],"kraken","standard")
    else:
        config["kraken_index_path"] = opj(config["resource_path"],"kraken","prebuilt", config["kraken_prebuilt"])
    if config["kraken_reduce_memory"]:
        config["kraken_params"]="--memory-mapping"
    else:
        config["kraken_params"]=""
