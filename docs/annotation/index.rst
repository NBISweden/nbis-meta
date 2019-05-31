Annotation of coding sequences
==============================

This section deals with how to annotate protein sequences predicted on
assembled contigs.

Protein family annotation
-------------------------
Annotation of predicted protein coding sequences is performed using the tools
`eggnog-mapper <https://github.com/jhcepas/eggnog-mapper>`_
`pfam_scan <https://www.ebi.ac.uk/Tools/pfa/pfamscan/>`_ (using the latest version of the PFAM database) and
`rgi <https://github.com/arpcard/rgi>`_ (Resistance Gene Identifier).

:code:`eggnog:` Set to True in order to annotate coding sequences with eggnog-mapper.
This also adds annotations for KEGG orthologs, modules and pathways.

:code:`pfam:` Set to True in order to add PFAM protein families to coding sequences.

:code:`rgi`: Set to True in order to add AMR (Antimicrobial Resistance) gene families to coding sequences.

:code:`rgi_params`: Determines how to run the rgi software. By default it's run in 'perfect' mode, using the diamond alignment tool.
See the `rgi GitHub pages <https://github.com/arpcard/rgi#rgi-main-usage-for-genomes-genome-assemblies-metagenomic-contigs-or-proteomes>_`
for more information.

.. note::
    The eggnog-mapper, pfam_scan and rgi programs all utilize separate conda environments. In order to use these tools you must run
    the workflow with the :code:`--use-conda` flag.

Taxonomic annotation
--------------------
Taxonomic annotation of protein sequences is performed by:
1. similarity searches against a reference database
2. parsing best hits for each protein query and linking hits to taxonomic IDs
3. assigning a taxonomic label to the lowest common taxa for the 'best hits'

The first step above is performed in this workflow using
`Diamond <https://github.com/bbuchfink/diamond/>`_ against a protein sequence
reference database. The workflow supports using the uniref50, uniref90, uniref100 or nr protein databases but you can
also use a preformatted custom database as long as it has been formatted with

:code:`taxonomic_annotation:` Set to True to add taxonomic annotations to contigs (and ORFs).


:code:`taxdb:` Can be one of 'nr', 'uniref50', 'uniref90', 'uniref100' or you can add a custom diamond database if you
format it with taxonomic ids (using the --taxonmap and --taxonnodes flags to diamond makedb)

:code:`diamond_threads`: How many threads to use for the diamond blast job.

:code:`taxonomy_min_len`: Minimum length of contigs to be included in the taxonomic annotation.

:code:`tango_params`: Additional parameters to use for the `tango <https://github.com/johnne/tango>`_ taxonomic assigner

:code:`taxonomy_ranks`: Taxonomic ranks to use in output.

.. note::
    If you want to use a custom preformatted protein database. Make sure to put it as 'diamond.dmnd' in the
    <resource_path> directory as specified in your configfile. Then run the workflow once with the :code:`-t` flag
    to update the timestamps on the existing files.

Generating required database files
==================================

Several steps of this workflow requires large database files that may take
a long time to download and format. If you want to can run the database
creation separately with the workflow (e.g. while you're waiting for real
data to arrive).

To create the databases needed for the protein annotation steps (eggnog, pfam, rgi and taxonomic) you can run::

    snakemake --configfile config.yaml --config eggnog=True pfam=True rgi=True taxonomic_annotation=True -np db

    <lots of text>

    Job counts:
        count   jobs
        1       db
        1       db_done
        1       download_eggnog
        1       download_pfam
        1       download_pfam_info
        1       download_uniref
        1       get_eggnog_version
        1       get_kegg_ortholog_info
        1       prepare_diamond_db_uniref
        1       prepare_taxfiles_uniref
        1       press_pfam
        11


.. note:: The 'db' target only includes database files used in preprocessing (SortMeRNA) and protein annotation.

.. seealso:: See the documentation for ways to create databases for `read-classification <http://nbis-metagenomic-workflow.readthedocs.io/en/latest/classification/index.html>`_.


