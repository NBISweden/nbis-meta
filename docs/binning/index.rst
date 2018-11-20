Binning metagenomic contigs
===========================

Assembled contigs can be grouped into so called 'genome bins' using information
about their nucleotide composition and their abundance profiles in several
samples. The rationale is that contigs from the same genome will have similar
nucleotide composition and will show up in similar abundances across
samples.

There are several programs that perform unsupervised clustering of contigs
based on composition and coverage. Currently, this workflow only uses
`MaxBin2 <https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html>`_
but other popular similar tools are `CONCOCT <https://github.com/BinPro/CONCOCT/>`_
and `MetaWatt <https://sourceforge.net/projects/metawatt/>`_.

The quality and phylogeny of bins can also be estimated as part of the workflow, using
`CheckM <https://github.com/Ecogenomics/CheckM>`_.

To perform the binning step of the workflow run the following::

    snakemake --use-conda --configfile <your-config-file> -p binning

.. seealso:: Check out the binning tutorial further down in this document.

Settings
--------
:code:`binning:` Set to True to perform binning on all assembly groups in the sample annotation file

Output from cross-mapping reads to contigs are saved in <results_path>/binning/map. Bam files are marked as
temporary in the workflow and thus the map_dir will not take up a huge amount of space once
the abundance files are prepared.

Output from the actual binning is placed in <results_path>/binning/bin. Each binned assembly will have a sub-folder under
the bin_dir directory. In each sub-folder a nucleotide fasta file will be generated for each genome
bin with the filename pattern "{assemblyGroup}.nnn.fasta" where 'nnn' is a number starting
at 001. The output folder will also contain a log file from the binning, a summary file
with some stats on each bin, a list of contigs that could not be binned as well as a list
of contigs that were too short (below the min_contig_length specified).

.. Note:: The completeness estimate in MaxBin2 does not work properly so don't believe the completeness values in the .summary file. Instead run the CheckM step of the workflow to get better estimates.

:code:`min_contig_length:` Minimum length of contigs to be binned. Shorter contigs will be filtered out.

:code:`binning_threads:` Number of threads to use for the MaxBin2 step.

:code:`checkm:` Set to True to analyze bins with Checkm.

Binning tutorial
----------------
When the binning step is included in the workflow all assembly groups
specified in the `sample annotation file <http://nbis-metagenomic-workflow.readthedocs.io/en/latest/configuration/sample_list.html>`_
will be binned and for each assembly the abundances of contigs will be
calculated by cross-mapping all samples in the sample annotation file.

This means that if you have many assemblies and/or samples there will be
a lot of jobs to run which can be a bit overwhelming. If you want you may
first run the `assembly <http://nbis-metagenomic-workflow.readthedocs.io/en/latest/assembly/index.html>`_
step and then bin individual assemblies.

Let's try this out using some supplied example data. The configuration
file :code:`examples/binning/config.yaml` is set up to use four mock communities created with this
`metagenomic-mocks <https://bitbucket.org/johnne/metagenomic-mocks>`_ repository.
The samples each contain 100,000 read-pairs sampled from the same 10
genomes but in varying proportions, roughly simulating the proportions
of the genomes in different body sites.

Have a look at the jobs that will be run in this example by doing a dry-run::

    snakemake --use-conda --configfile examples/binning/config.yaml -np

There are a total of 204 jobs to run here. Results will be placed under :code:`results/examples/binning`.

To also analyze the quality of generated bins using CheckM you can either set :code:`checkm:True` in your config file
or from the command line::

    snakemake --use-conda -j 4 --configfile examples/binning/config.yaml --config checkm=True -np

.. Important:: Notice the :code:`--use-conda` flag to snakemake. This is required when running checkm as that software uses a different conda environment.

To run the binning example using 4 cores with bin QC using CheckM simply do::

    snakemake --use-conda -j 4 --configfile examples/binning/config.yaml --config checkm=True -p

Splitting up the example
^^^^^^^^^^^^^^^^^^^^^^^^

If you want, you can first generate the assemblies which will be used for
binning. Simply run::

    snakemake --use-conda -j 4 --configfile examples/binning/config.yaml -p assembly

Then bin each assembly separately::

    snakemake --use-conda --configfile examples/binning/config.yaml -p results/examples/binning/bins/$ASSEMBLY/$ASSEMBLY.summary

Here :code:`$ASSEMBLY` should be substituted for each of the assemblies generated.

Produce CheckM quality plot::

    snakemake --use-conda --configfile examples/binning/config.yaml -np results/examples/checkm/bin_qa_plot.png
