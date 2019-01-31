Binning metagenomic contigs
===========================

Assembled contigs can be grouped into so called 'genome bins' using information
about their nucleotide composition and their abundance profiles in several
samples. The rationale is that contigs from the same genome will have similar
nucleotide composition and will show up in similar abundances across
samples.

There are several programs that perform unsupervised clustering of contigs
based on composition and coverage. Currently, this workflow uses
`MaxBin2 <https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html>`_ and/or
`CONCOCT <https://github.com/BinPro/CONCOCT/>`_. The quality and phylogeny of bins can be estimated using
`CheckM <https://github.com/Ecogenomics/CheckM>`_.

.. Note:: The minimum required version of CONCOCT is 1.0.0 which only works on Linux so far. However, if you manage to compile CONCOCT on OSX you can use it with the workflow by installing it into the nbis-meta environment path and then omit the :code:`--use-conda` flag to snakemake.


To perform the binning step of the workflow run the following::

    snakemake --use-conda --configfile <your-config-file> -p binning

.. seealso:: Check out the binning tutorial further down in this document.

Settings
--------
:code:`maxbin:` Set to True to use MaxBin2.

:code:`concoct:` Set to True to use CONCOCT.

:code:`min_contig_length:` Minimum contig length to include in binning. Shorter contigs will be ignored. By setting
multiple values here you can make the workflow run the binning steps multiple times with different minimum lengths.
Output is stored under :code:`results/maxbin/<min_contig_length>/` and :code:`results/concoct/<min_contig_length>/`.

:code:`maxbin_threads:` Number of threads to use for the MaxBin2 step.

:code:`concoct_threads:` Number of threads to use for the CONCOCT step.

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

Results will be placed under :code:`results/examples/`.

Splitting up the example
^^^^^^^^^^^^^^^^^^^^^^^^

If you want, you can first generate the assemblies which will be used for
binning. Simply run::

    snakemake -j 4 --configfile examples/binning/config.yaml -p assembly

Then run the binning steps::

    snakemake --use-conda -j 4 --configfile examples/binning/config.yaml -p binning
