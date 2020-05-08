# NBIS-Metagenomics
A workflow for metagenomic projects

[![Python 3.7.6](https://img.shields.io/badge/python-3.7.6-brightgreen.svg)](https://www.python.org/downloads/release/python-376/)
[![Snakemake 5.11.2](https://img.shields.io/badge/snakemake-5.11.2-brightgreen.svg)](https://img.shields.io/badge/snakemake-5.11.2)
![CI](https://github.com/NBISweden/nbis-meta/workflows/CI/badge.svg?branch=master)
![Docker](https://img.shields.io/docker/pulls/nbisweden/nbis-meta)


## Overview
A [snakemake](http://snakemake.readthedocs.io/en/stable/) workflow for
paired- and/or single-end metagenomic data.

You can use this workflow for _e.g._:

- **read-trimming and QC**
- **taxonomic classification**
- **assembly**
- **functional and taxonomic annotation**
- **metagenomic binning**

See the [Wiki-pages](https://github.com/NBISweden/nbis-meta/wiki) for 
instructions on how to run the workflow.

## Installation

### From GitHub
Checkout the latest version of this repository:

```
git clone https://github.com/NBISweden/nbis-meta.git
```

This workflow runs on a minimal conda environment with `snakemake`, `biopython`
and `pandas`. To get started use the [conda](http://anaconda.org) environment
file `environment.yml`. 

(If you don't have conda installed already refer to the
[conda installation instructions](https://conda.io/docs/user-guide/install/index.html).)

To install the workflow environment, run:

```
conda env create -f environment.yml
```

Then activate the environment with:

```
conda activate nbis-meta
```

### With Docker

To pull the latest Docker image with all dependencies and source code from
DockerHub, run:

```bash
docker pull nbisweden/nbis-meta
```

See the Wiki for instructions on how to run the Workflow with Docker.