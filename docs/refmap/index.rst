Reference based mapping
=======================

This part of the workflow performs mapping of preprocessed reads against
a set of reference sequences and reports on the abundance of genomes
and their coding sequences in your samples.

If you're working with samples that you are fairly certain contain
organisms with representative genome references this type of analysis
may suit you.

Prefiltering
------------
The workflow will first use the read-based classifier `Centrifuge <https://github.com/infphilo/centrifuge>`_
to make a rough estimate of the genomes present in your samples. The
main advantage of this classifier is that the database doesn't require
a lot of storage (the full refseq database of bacteria, archaea and viruses
only takes up ~5 GB) so it's good for doing a quick filtering run.

Take a look at `How to configure the Centrifuge database <http://nbis-metagenomic-workflow.readthedocs.io/en/latest/classification/index.html#id1>`_)
for more details.

The prefiltering step classifies all samples against the centrifuge database
and filters out taxa that have at least :code:`refmap_min_read_count:` mapped. The
default is 50 reads but this can be changed in the config file or on the command line by adding
:code:`--config refmap_min_read_count <min_reads>`.

Genomic data for genomes meeting the threshold in any of the samples is
downloaded and a bowtie2 index is built.

Mapping
-------
The workflow proceeds with mapping of reads against the set of filtered
genomes using bowtie2. Read mappings are then counted for each coding sequence
defined in the reference GFF files that has a 'protein_id' in the attributes field.

Output
------
Raw counts and normalized abundances (TPM) of the coding sequences are reported
for all samples as well as a file of total genome coverage for the filtered genomes.

PCA plots of samples based on raw and normalized coding sequence abundance as well as
genome coverage are also generated.

Example
-------
To try out the reference based mapping you can use the example data under
:code:`examples/data/` and the corresponding configuration file :code:`examples/refmap/config.yaml`.

To speed up the building of the centrifuge database this example limits the
organisms to include to a set of 28 taxonomy ids specified in :code:`examples/refmap/taxids`.

Simply run::

    snakemake --configfile examples/refmap/config.yaml -p -j 4

to run the workflow with 4 cores (modify the :code:`-j` parameter to change this).