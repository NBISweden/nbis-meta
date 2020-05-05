# NBIS-Metagenomics
A workflow for metagenomic projects

[![Python 3.7.6](https://img.shields.io/badge/python-3.7.6-blue.svg)](https://www.python.org/downloads/release/python-376/)
[![Snakemake 5.11.2](https://img.shields.io/badge/snakemake-5.11.2-brightgreen.svg?style=flat-square)](https://img.shields.io/badge/snakemake-5.11.2)
![Docker](https://img.shields.io/docker/pulls/nbisweden/nbis-meta)

## Overview
A [snakemake](http://snakemake.readthedocs.io/en/stable/) workflow for
paired- and/or single-end metagenomic data.

You can use this workflow for, _e.g._:

- **read-trimming**
- **taxonomic classification**
- **assembly**
- **functional and taxonomic annotation**
- **metagenomic binning**

## Installation

### From GitHub
Checkout the latest version of this repository:

```
git clone https://github.com/NBISweden/nbis-meta.git
```

Change directory:

```
cd nbis-meta
```

#### Install the required software
All the software needed to run this workflow is included as a
[conda](http://anaconda.org) environment file. You will need to
[install conda](https://conda.io/docs/user-guide/install/index.html)
before installing the required software for this workflow.

To create the environment `nbis-meta` use the supplied
[environment.yaml](envs/environment.yaml) file found in the *envs/*
folder.

First create the environment using the supplied file:

```
conda env create -f envs/environment.yml
```

Activate the environment using:

```
conda activate nbis-meta
```

#### Configure workflow for the SLURM Workload Manager (e.g. Uppmax)
If you are going to run the workflow on a compute cluster such as
[Uppmax](https://uppmax.uu.se/) you can make use of the snakemake SLURM
profile created by [Per Unneberg](https://github.com/percyfal).

**Create a directory to hold the profile:**

```
mkdir profiles
```

**Install and activate cookiecutter with conda:**

```
mkdir envs/cookiecutter
conda env create -f envs/cookiecutter.yml -p envs/cookiecutter
conda activate envs/cookiecutter
```

**Install the snakemake profile into the profiles/ directory:**

```
cookiecutter -o profiles https://github.com/Snakemake-Profiles/slurm.git
```

You will be prompted for the account to charge compute hours to as well
as the default partition. Once this is done you can run snakemake on
the cluster using `snakemake --profile profiles/slurm -j 100`.

## Documentation
See the [documentation](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/index.html) for instructions on how to run the pipeline.
