==========================
NBIS Metagenomics Workflow
==========================
A snakemake workflow for preprocessing and analysis of metagenomic
datasets. Handles single- and/or paired-end data and runs on a local UNIX/LINUX
laptop, or in a cluster environment.

The source code is available at
`GitHub <https://github.com/NBISweden/nbis-meta>`_ and is being developed as
part of the `NBIS <http://nbis.se>`_ bioinformatics infrastructure.


Quick start
===========

From GitHub
-----------

Clone the repository::

    git clone git@github.com:NBISweden/nbis-meta.git
    cd nbis-meta

Install and activate the conda environment::

    conda env create -f environment.yml
    conda activate nbis-meta

**Done!**

To run the workflow on some example data simply type::

    snakemake --use-conda -j 4 -p


With Docker
-----------

Pull the latest Docker image containing everything needed to run the workflow::

    docker pull nbisweden/nbis-meta

Run the workflow in a container (mounting the current directory)::

    docker run -v $(pwd):/analysis nbisweden/nbis-meta

.. toctree::
    :maxdepth: 1
    :caption: Getting started
    :hidden:

    starting

.. toctree::
    :maxdepth: 1
    :caption: Preprocessing
    :hidden:

    preprocessing

.. toctree::
    :maxdepth: 1
    :caption: Read-based classification
    :hidden:

    read-based

.. toctree::
    :maxdepth: 1
    :caption: De-novo assembly
    :hidden:

    denovo