NBIS Metagenomics Workflow
====================================

.. toctree::
    :maxdepth: 1
    :caption: Setup

    configuration/index

.. toctree::
    :maxdepth: 1
    :caption: Preprocessing

    preprocessing/index

.. toctree::
    :maxdepth: 1
    :caption: Read-based

    classification/index
    profiling/index
    refmap/index

.. toctree::
    :maxdepth: 1
    :caption: De-novo assembly

    assembly/index
    annotation/index
    binning/index

Overview
--------
This is a snakemake workflow for rocessing and analysis of metagenomic
datasets. It can handle single- and paired-end data and can run on a
local laptop with either Linux or OSX, or in a cluster environment.

The source code is available at `BitBucket <https://bitbucket.org/scilifelab-lts/nbis-meta>`_ and
is being developed as part of the `NBIS <http://nbis.se>`_ bioinformatics infrastructure.

Installation
------------
**1. Clone the repository**
Checkout the latest version of this repository (to your current directory)::

    git clone https://bitbucket.org/scilifelab-lts/nbis-meta

**2. Install the required software**
All the software needed to run this workflow is included as a
`Conda <http://anaconda.org>`_ environment file. To create the
environment :code:`nbis-meta` use the supplied :code:`envs/environment.yaml` file.

First create the environment using the supplied file::

    mkdir envs/nbis-meta
    conda env create -f envs/environment.yaml -p envs/nbis-meta

Next, add this directory to the envs_dirs in your conda config (this is to simplify
activation of the environment and so that the full path of the
environment installation isn't shown in your bash prompt)::

    conda config --add envs_dirs $(pwd)/envs/

Activate the environment using::

    conda activate envs/nbis-meta

**You are now ready to start using the workflow!**

.. note:: If you plan on using the workflow in a cluster environment running the SLURM workload manager (such as Uppmax) you should configure the workflow with the SLURM snakemake profile. `See the documentation <https://nbis-metagenomic-workflow.readthedocs.io/en/latest/configuration/index.html#how-to-run-on-uppmax-hebbe-snic-resources>`_.
