# NBIS-Metagenomics
A workflow for metagenomic projects

[![Documentation Status](https://readthedocs.org/projects/nbis-metagenomic-workflow/badge/?version=latest)](http://nbis-metagenomic-workflow.readthedocs.io/en/latest/?badge=latest)
[![Python 3.5](https://img.shields.io/badge/python-3.5-blue.svg)](https://www.python.org/downloads/release/python-350/)
[![Snakemake 4.4](https://img.shields.io/badge/snakemake-%E2%89%A54.4.0-brightgreen.svg?style=flat-square)](https://img.shields.io/badge/snakemake-%E2%89%A54.4.0-brightgreen.svg?style=flat-square)

## Overview
This is a [snakemake](http://snakemake.readthedocs.io/en/stable/) workflow that processes paired-end and/or single-end metagenomic samples.

Potential analyses include:

- **read-based taxonomic classification**
- **assembly**
- **functional and taxonomic annotation** of coding sequences
- **genome binning** of assembled contigs

## Installation

### Clone the repository
Checkout the latest version of this repository:

```
git clone https://github.com/NBISweden/nbis-meta.git
```

Change directory:

```
cd nbis-meta
```

### Install the required software
All the software needed to run this workflow is included as a
[conda](http://anaconda.org) environment file. You will need to
[install conda](https://conda.io/docs/user-guide/install/index.html)
before installing the required software for this workflow.

To create the environment `nbis-meta` use the supplied
[environment.yaml](envs/environment.yaml) file found in the *envs/*
folder.

First create the environment using the supplied file:

```
mkdir -p envs/nbis-meta
conda env create -f envs/environment.yaml -p envs/nbis-meta
```

This creates the `nbis-meta` environment inside the `envs/` directory and
installs the environment there.

Next, add this directory to the `envs_dirs` in your conda config (this is to 
simplify activation of the environment).

```
conda config --add envs_dirs $(pwd)/envs/
```

Activate the environment using:

```
conda activate nbis-meta
```

### Configure workflow for the SLURM Workload Manager (e.g. Uppmax)
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
conda env create -f envs/cookiecutter.yaml -p envs/cookiecutter
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
