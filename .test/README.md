# Testing

This directory contains config files, resources and scripts used to test the
workflow.

**config/**

Files under `config/` are used by the different test steps of the github
actions testing.

**data/**

The fasta file at `data/uniref100.fasta` was created by:

1. searching [uniprot](https://uniprot.org) for the 5 taxids used to generate
the [synthetic metagenome](https://zenodo.org/record/3737112#.XsUQncZ8LOQ) that
this workflow uses for testing. This resulted in 8,122 identified sequences.
2. mapping the sequences to their UniRef100 id via the Retrieve/ID mapping tool
at uniprot, followed by downloading of the reference sequences.
3. subsampling 100 sequences from the downloaded fastafile using `seqtk

During testing the fasta file is used to build and query a diamond database
using `tango`.
