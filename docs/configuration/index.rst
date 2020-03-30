Installation
============
**1. Clone the repository**
Checkout the latest version of this repository::

    git clone https://bitbucket.org/scilifelab-lts/nbis-meta

Change directory::

    cd nbis-meta

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

Quick start guide
=================

Using the examples
-------------------
The workflow comes with a set of example datasets found in
the :code:`examples/data` folder. By default the workflow uses the
configuration specified in the :code:`config.yaml`
file to do **preprocessing** and **assembly** of the example data.

To see what the workflow is configured to do in the default setup run::

    snakemake -np

This performs a 'dry run' (specified by the :code:`-n` flag) where nothing is
actually executed just printed to the terminal. Try it out by removing
the :code:`-n`.

To familiarize yourself with ways to run the workflow you can make use
of the example directories found under :code:`examples/`.

For instance, to run the centrifuge example (Read classification and
taxonomic profiling using Centrifuge) do::

    snakemake --configfile examples/centrifuge/config.yaml

See more about the examples here: `binning`_, `centrifuge`_

.. _binning: https://nbis-metagenomic-workflow.readthedocs.io/en/latest/binning/index.html#binning-tutorial
.. _centrifuge: http://nbis-metagenomic-workflow.readthedocs.io/en/latest/classification/index.html#example-run-with-centrifuge

Using your own data
-------------------
1. Create your `sample_list <http://nbis-metagenomic-workflow.readthedocs.io/en/latest/configuration/sample_list.html>`_
2. Configure the workflow to your needs. Make sure to updated the :code:`sample_list:` path in the config.

.. note:: We recommend that you copy the config.yaml file to some other file, e.g. 'myconfig.yaml', and make your changes in the copy.

You are now ready to start the workflow using::

    snakemake --configfile myconfig.yaml

You can have multiple configuration files (e.g. if you want to run the
workflow on different sets of samples or with different parameters).


General configuration settings
==============================

Here you'll find information on how to configure the pipeline.

The file :code:`config.yaml` contains all the parameters for the pipeline.
The most general settings are explained here. For more specific settings
regarding the different workflow steps you may refer to the corresponding documentation.

.. note:: The only required config setting is :code:`sample_list:`. Without this, the workflow cannot run but essentially your config file may contain only this parameter.

:code:`workdir:` Working directory for the snakemake run. All paths are evaluated relative to this directory.

:code:`sample_list:` should point to a file with information about samples. See the documentation on :ref:`sample-list` for more details on the format of this file.

:code:`taxdb`: Path to where taxonomic database files from NCBI are stored.

:code:`results_path:` this directory will contain reports, figures, count tables and other types of
aggregated and processed data. Output from the different steps will be created in
sub-directories. For instance, assemblies will be in a sub-directory called 'assembly'

:code:`intermediate_path:` this directory contains files which might have required a substantial computational
effort to generate and which might be useful for further analysis. They are not deleted
by Snakemake once used by downstream rules. This could be for example bam files.

:code:`temp_path:` this directory contains files which are either easily regenerated or redundant. This
could be for example sam files which are later compressed to bam.

:code:`scratch_path:` Similar to temp_path but primarily used on cluster resources (e.g. Uppmax)
in order to speed up writes to disk by utilizing local storage for nodes.

:code:`resource_path:` Path to store database files needed for the workflow.

:code:`report_path:` Path to store report files such as assembly and mapping stats and plots

.. _sample-list:

The sample list file
====================
The :code:`sample_list:` config parameter should point to a file containing information about the samples and respective
files to be used as input for the pipeline.

The sample list file specifies where the input data for the workflow is located on your computer or
the compute cluster you are running on.

Below is an example of what the sample file may look like. This
example table is also found in :code:`samples/sample_annotation_example.tab`.
This example uses data from the repository `Metagenomic Mocks <https://bitbucket.org/johnne/metagenomic-mocks>`_. Each sample contains 100 000 read-pairs.

===============  =====  ===================  =================================================   =================================================
sampleID         runID  assemblyGroup            fileName                                           pair
===============  =====  ===================  =================================================   =================================================
 anterior_nares    1     anterior_nares,all   examples/data/anterior_nares_100000_R1.fastq.gz     examples/data/anterior_nares_100000_R2.fastq.gz
 buccal_mucosa     1     buccal_mucosa,all    examples/data/buccal_mucosa_100000_R1.fastq.gz      examples/data/buccal_mucosa_100000_R2.fastq.gz
 retr_crease       1     retr_crease,all      examples/data/retr_crease_100000_R1.fastq.gz        examples/data/retr_crease_100000_R2.fastq.gz
 stool             1     stool,all            examples/data/stool_100000_R1.fastq.gz              examples/data/stool_100000_R2.fastq.gz
===============  =====  ===================  =================================================   =================================================

**The sampleID and runID columns:**

The *sampleID* column allows you to designate a sample ID for each set of sequences while the *runID* column can
be used to designate e.g. technical replicates of samples. These two columns together form a unique tag for each
sequence set. If there is only one sequencing run per sampleID you may leave the runID column empty or simply
fill in a '1'.

**The assemblyGroup column:**

The *assemblyGroup* column allows you to group together samples (and/or individual
sample runs) into assembly groups. A single sample/run combination can be grouped into multiple assembly groups by
specifying comma separated assembly group names in this field. In the example above each sample has been assigned to
an individual assembly as well as a co-assembly named 'all' which will contain all samples. Running the workflow
with this file will produce five assemblies in total (named 'anterior_nares', 'buccal_mucosa', 'retr_crease', 'stool'
and 'all).

**The fileName and pair columns:**

These two columns specify file paths for sequences in the (gzipped) fastq format.
For paired end data the *fileName* column points to *forward* read file and the *pair* column points to the
corresponding *reverse* read file. For single end data only the *fileName* column is
used.


How to run on UppMax/Hebbe (SNIC resources)
===========================================
The recommended way to run this workflow on a SLURM cluster such as Uppmax is to install the
`SLURM profile <https://github.com/Snakemake-Profiles/slurm>`_ for snakemake.

To do so you will need cookiecutter which you can install using the supplied environment file::

    mkdir envs/cookiecutter
    conda env create -f envs/cookiecutter.yaml -p envs/cookiecutter

Then activate the cookiecutter environment and deploy the profile::

    conda activate envs/cookiecutter
    mkdir profiles
    cookiecutter -o profiles https://github.com/Snakemake-Profiles/slurm.git

You will be prompted to add some information such as account number, partition etc. You can leave some of these fields
blank but should at least fill out the account number (e.g. snic2017-1-234 on SNIC resources).
Below is a recommended example::

    account []: snic2017-1-234 # Use your actual account number!
    error []:
    output []:
    partition []: core
    profile_name [slurm]: slurm
    Select submit_script:
    1 - slurm-submit.py
    2 - slurm-submit-advanced.py
    Choose from 1, 2 (1,2) [1]: 2

You can now run the workflow in the cluster environment using::

    snakemake --profile profiles/slurm -j 100 -np

The :code:`-j 100` flag tells snakemake to have at most 100 jobs submitted to the SLURM queue at the same time.
