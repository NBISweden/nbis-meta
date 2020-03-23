==========================
NBIS Metagenomics Workflow
==========================

.. toctree::
    :maxdepth: 1
    :caption: Setup
    :hidden:

    configuration/index

.. toctree::
    :maxdepth: 1
    :caption: Preprocessing
    :hidden:

    preprocessing/index

.. toctree::
    :maxdepth: 1
    :caption: Read-based classification
    :hidden:

    classification/index

.. toctree::
    :maxdepth: 1
    :caption: De-novo assembly
    :hidden:

    assembly/index
    annotation/index
    binning/index

========
Overview
========
This is a snakemake workflow for preprocessing and analysis of metagenomic
datasets. It can handle single- and paired-end data and can run on a
local laptop with either Linux or OSX, or in a cluster environment.

The source code is available at
`GitHub <https://github.com/NBISweden/nbis-meta>`_ and is being developed as
part of the `NBIS <http://nbis.se>`_ bioinformatics infrastructure.

============
Installation
============

From GitHub
-----------

**1. Clone the repository**

Checkout the latest version of this repository (to your current directory)::

    git clone git@github.com:NBISweden/nbis-meta.git
    cd nbis-meta

Change directory::

    cd nbis-meta

**2. Install the required software**

The workflow is structured around a default core set of rules and several
additional rules that can be set my modifying the config settings. Because of
this modular structure the software requirements are divided in several conda
files. This ensures that the size of the software environment is only as big
as it needs to be for your analysis. However, it also means that snakemake
**has to be run with the** :code:`--use-conda` **flag**.

The bare-minimum requirements for starting the workflow can be installed with
the :code:`environment.yml` file::

    conda env create -f environment.yml


Activate the environment using::

    conda activate envs/nbis-meta

**You are now ready to start using the workflow!**

To run the workflow on some example data simply type::

    snakemake --use-conda -j 4 -p

From Docker
-----------

You can also pull the latest Docker image which has all the code and
requirements needed to run the workflow::

    docker pull nbisweden/nbis-meta:latest

To run a container with the image and get an interactive shell, run::

    docker run --rm -it -v $(pwd):/analysis nbisweden/nbis-meta:latest /bin/bash

.. note::
    If you plan on using the workflow in a cluster environment running the
    SLURM workload manager (such as Uppmax) you should configure the workflow
    with the SLURM snakemake profile.
    `See the documentation <https://nbis-metagenomic-workflow.readthedocs.io/en/latest/configuration/index.html#how-to-run-on-uppmax-hebbe-snic-resources>`_.
