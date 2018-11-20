# Documentation

Contents
--------
* [Annotation](annotation/index.md)
* [Assembly](assembly/index.md)
* [Binning](binning/index.md)
* [Classifying reads](classification/index.md)
* [Configuration](configuration/index.md)
* [Preprocessing](preprocessing/index.md)
* [Reference-based mapping](refmap/index.md)

## Overview
This is a snakemake workflow for rocessing and analysis of metagenomic
datasets. It can handle single- and paired-end data and can run on a
local laptop with either Linux or OSX, or in a cluster environment.

The source code is available at [BitBucket](https://bitbucket.org/scilifelab-lts/lts-workflows-sm-metagenomics)
 and is being developed as part of the [NBIS](http://nbis.se)
 bioinformatics infrastructure.

```eval_rst
.. image:: http://nbis.se/assets/img/logos/nbislogo-green-txt.svg
    :width: 200
    :alt: NBIS logo
```

## Installation

### 1. Clone the repository
Checkout the latest version of this repository (to your current directory):

```
git clone https://bitbucket.org/scilifelab-lts/lts-workflows-sm-metagenomics
```

### 2. Install the required software
All the software needed to run this workflow is included as a
[conda](http://anaconda.org) environment file. To create the
environment `sm-meta` use the supplied
[environment.yaml](lts_workflows_sm_metagenomics/envs/environment.yaml) file
found in the *lts_workflows_sm_metagenomics/envs/* folder.

First change directory to the lts_workflows_sm_metagenomics subdir:

```
cd lts-workflows-sm-metagenomics/lts_workflows_sm_metagenomics
```

Then create the environment using the supplied file:

```
conda env create -f envs/environment.yaml
```

Activate the environment using:

```
source activate sm-meta
```

**You are now ready to start using the workflow!**

## Quick start guide
The workflow comes with a set of example datasets found in
the `examples/data` folder. By default the workflow uses the
configuration specified in the `default_configuration_parameters.yaml`
file to do **preprocessing** and **assembly** of the example data.

To see what the workflow is configured to do in the default setup run:

```
snakemake -np
```

This performs a 'dry run' (specified by the `-n` flag) where nothing is
actually executed just printed to the terminal. Try it out by removing
the `-n`.

### Using the examples
To familiarize yourself with ways to run the workflow you can make use
of the example directories found under `examples/`.

For instance, to run the centrifuge example (Read classification and
taxonomic profiling using Centrifuge) do:

```
snakemake --configfile examples/centrifuge/config.yaml centrifuge_classify
```

See more about the examples here: [binning](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/binning/index.html#binning-tutorial),
 [centrifuge](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/classification/index.html#example-run-with-centrifuge),
 [reference-based mapping](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/refmap/index.html#example)

### Using your own data
1. Create your [sample_list](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/configuration/sample_list.html).
2. Configure the workflow to your needs. Make sure to updated the
  `sample_list:` path in the config. *We recommend that you copy
the default_configuration_parameters.yaml to some other file, e.g.
'myconfig.yaml', and make your changes in the copy.*

You are now ready to start the workflow using:

```
snakemake --configfile myconfig.yaml
```

You can have multiple configuration files (e.g. if you want to run the
workflow on different sets of samples or with different parameters).

## Pipeline steps
If you don't want to run all the steps of this pipeline in one go you can
specify [targets](http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#targets)
which will cause snakemake to run only certain jobs. The available targets are
listed below. Click on the links to see a more detailed explanation and
related settings.

* __[annotation](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/annotation/index.html)__: Annotation of coding sequences found on assembled contigs.
* __[assembly](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/assembly/assembly.html)__: De-novo assembly of reads into contigs.
* __[binning](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/binning/binning.html)__: Binning of metagenomic contigs.
* __[centrifuge_classify](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/classification/index.html)__: Read classification and taxonomic profiling using Centrifuge.
* __[db](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/annotation/db.html)__: Setup of databases for protein sequence annotation
* __[kraken_classify](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/classification/index.html)__: Read classification and taxonomic profiling using Kraken.
* __[preprocess](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/preprocessing/preprocess.html)__: Preprocessing of reads.
* __[refmap](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/refmap/reference_based_mapping.html)__: Reference-based read mapping.