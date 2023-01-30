# NBIS-Metagenomics
A workflow for metagenomic projects

[![Python 3.9](https://img.shields.io/badge/python-3.9-brightgreen.svg)](https://www.python.org/downloads/release/python-390/)
[![Snakemake 7.0.0](https://img.shields.io/badge/snakemake-7.0.0-brightgreen.svg)](https://img.shields.io/badge/snakemake-7.0.0)
![CI](https://github.com/NBISweden/nbis-meta/workflows/CI/badge.svg?branch=main)
![Docker](https://img.shields.io/docker/pulls/nbisweden/nbis-meta)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

## Overview
A [snakemake](http://snakemake.readthedocs.io/en/stable/) workflow for
paired- and/or single-end metagenomic data.

You can use this workflow for _e.g._:

- **read-trimming and QC**
- **taxonomic classification**
- **assembly**
- **functional and taxonomic annotation**
- **metagenomic binning**

See the [documentation](https://nbis-metagenomic-workflow.readthedocs.io/en/latest/#) for
instructions on how to run the workflow.

## Installation

### From GitHub
1. Checkout the latest version:

```
git clone https://github.com/NBISweden/nbis-meta.git
```

or download a tarball of the latest release from the [release page](https://github.com/NBISweden/nbis-meta/releases).

2. Install and activate the workflow environment:

```
mamba env create -f environment.yml
conda activate nbis-meta
```

### From DockerHub

To pull the latest Docker image with all dependencies and source code from
DockerHub, run:

```bash
docker pull nbisweden/nbis-meta
```

See the Wiki for instructions on how to run the Workflow with Docker.
