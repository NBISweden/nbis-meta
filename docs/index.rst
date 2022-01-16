.. toctree::
    :hidden:

    setup
    howto
    uppmax

########
Overview
########

nbis-meta is a snakemake workflow for metagenomics projects

**********
Installing
**********

From GitHub
===========

To start using the workflow either clone the latest version of the repo, run:

.. code-block:: bash

    git clone https://github.com/NBISweden/nbis-meta.git

or download the latest release from the `release page <https://github.com/NBISweden/nbis-meta/releases>`_
and extract the archive.

Then change directory into the ``nbis-meta`` folder and create the core conda
environment:

.. code-block:: bash

    cd nbis-meta
    conda env create -f environment.yml

.. note:: mamba instead of conda

    ``mamba`` is a faster replacement for conda. Give it a try by installing it from
    the conda-forge channel: ``conda install -c conda-forge mamba``.
    You can then run ``mamba env create -f environment.yml``.

From DockerHub
==============

You may also pull the latest Docker image:

.. code-block:: bash

    docker pull nbisweden/nbis-meta:latest


************
What's next?
************

You are now ready to start using the workflow!

* for information on how to prepare necessary files see :doc:`Getting-started <setup>`
* then check out the :doc:`How-to page <howto>` for more info on how to run the workflow

*****************
Workflow overview
*****************

Preprocessing
=============

This workflow can perform preprocessing of paired- and/or single-end whole-genome
shotgun metagenomic data (in fastq-format) using *e.g.*:

* Trimmomatic (adapter/quality trimming)
* Cutadapt (adapter trimming)
* SortMeRNA (rRNA filtering)
* Fastuniq (de-duplication)
* FastQC and MultiQC (read QC and report generation)

Downstream analysis
===================

Read-based classification
-------------------------

Preprocessed reads can be used for taxonomic classification and profiling using tools such as:

* Kraken2
* Centrifuge
* MetaPhlAn3

producing taxonomic profiles of the samples, as well as interactive krona plots.

Assembly-based analysis
-----------------------

Preprocessed reads can also be assembled and analyzed further using tools such as:

* Megahit/MetaSPADES (for metagenomic assembly)
* prodigal (gene calling)
* pfam_scan, eggnog-mapper, Resistance Gene Identifier (protein level annotations)
* bowtie2 (mapping reads to contigs)
* featureCounts (assigning and counting mapped reads)
* edgeR + metagenomeSEQ (normalization of read counts for genes/features)
* contigtax + sourmash (taxonomic assignments)
* metabat2, CONCOCT, MaxBin2 (metagenomic binning)
* checkm (genome bin QC)
* GTDB-TK (genome bin phylogenetic assignments)
* fastANI (genome bin clustering)
