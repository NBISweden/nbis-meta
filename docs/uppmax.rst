##############################
Running the workflow on Uppmax
##############################

If you will be running the workflow on the `Uppmax <https://uppmax.uu.se/>`_
HPC clusters here are some helpful tips.

.. toctree::

**************
Uppmax profile
**************

The workflow comes with a pre-configured profile for the Uppmax cluster. To use
it, first make sure you set your SLURM account to the configuration file with
the ``slurm_acount`` parameter in your config, then run the workflow with
``--profile uppmax`` in the command line call, *e.g.*:

.. code-block::bash

    snakemake --profile uppmax --configfile myconfig.yaml

*************************
Setting up database files
*************************

All resources required to run this workflow are automatically downloaded as
needed when the workflow is executed. However, some of these (often large) files
are already installed in a central location on the system at `/sw/data`. This
means that you can make use of them for your workflow runs thereby saving you
time and reducing overall disk-usage for the system.

First create a ``resource/`` directory inside the directory where you will be
running the workflow::

    mkdir resources

eggNOG database
===============

The ``5.0`` version of the eggNOG database is installed in a central location on
Uppmax under ``/sw/data/uppnex/eggNOG/5.0``. To make use of this database
simply symlink the files in the directory into a directory ``resources/eggnog-mapper``:

.. code-block:: bash

    mkdir resources/eggnog-mapper
    ln -s /sw/data/eggNOG/5.0/$CLUSTER/eggnog{.db,_proteins.dmnd} resources/eggnog-mapper/
    head -1 /sw/data/eggNOG/eggNOG_data-5.0_install-README.md > resources/eggnog-mapper/eggnog.version
    touch resources/eggnog-mapper/download.log

Pfam database
=============

To use the Pfam database from the central location, create a ``pfam``
sub-directory under ``resources`` and link the necessary files from the central
location, run the following:

.. code-block:: bash

    mkdir resources/pfam
    version="35.0"
    ln -s /sw/data/Pfam/$version/Pfam-A.hmm* resources/pfam/
    cat /sw/data/Pfam/$version/Pfam.version > resources/pfam/Pfam-A.version

This installs the necessary files for release ``31.0``. Check the directories
under ``/sw/data/Pfam/`` to see available release versions.

Kraken database
===============

For ``kraken`` there are a number of databases installed under ``/sw/data/Kraken2``.
Snapshots of the ``standard``, ``nt``, ``rdp``, ``silva`` and ``greengenes``
indices are installed on a monthly basis. To use the latest version of the
standard index, do the following:

1. Create a sub-directory and link the index files

.. code-block:: bash

    mkdir -p resources/kraken/standard
    ln -s /sw/data/Kraken2_data/latest/*.k2d resources/kraken/standard/

2. From a reproducibility perspective it's essential to keep track of when the
index was created, so generate a version file inside your kraken directory by
running:

.. code-block:: bash

    file /sw/data/Kraken2_data/latest | egrep -o "[0-9]{8}\-[0-9]{6}" > resources/kraken/standard/kraken.version

3. Modify your config file so that it contains:

.. code-block:: yaml

    classification:
        kraken: True

    kraken:
        # generate the standard kraken database?
        standard_db: True

Non-standard databases
----------------------

Using any of the other non-standard databases from the central location is also
a simple process, *e.g.* for the latest SILVA index:

.. code-block:: bash

    mkdir -p resources/kraken/silva
    ln -s /sw/data/Kraken2_data/latest_silva/*.k2d resources/kraken/silva/
    file /sw/data/Kraken2_data/latest_silva | egrep -o "[0-9]{8}\-[0-9]{6}" > resources/kraken/silva/kraken.version

Then update your config file with:

.. code-block:: yaml

    kraken:
      standard_db: False
      custom: "resources/kraken/silva"

Centrifuge database
===================

There are several centrifuge indices on Uppmax located at ``/sw/data/uppnex/Centrifuge-indices/20180720/``,
but keep in mind that they are from 2018.

The available indices are: ``p_compressed``, ``p_compressed+h+v``, ``p+h+v``,
and ``p+v`` (see the `centrifuge manual <https://ccb.jhu.edu/software/centrifuge/manual.shtml>`_
for information on what these contain).

To use the ``p_compressed`` index, run:

.. code-block:: bash

    mkdir -p resources/centrifuge/p_compressed/
    ln -s /sw/data/Centrifuge-indices/20180720/p_compressed.*.cf resources/centrifuge/p_compressed/

Then update your config file to contain:

.. code-block:: yaml

    classification:
        centrifuge: True

    centrifuge:
        custom: "resources/centrifuge/p_compressed/p_compressed"

GTDB
===================

To use the centrally installed Genome Taxonomy Database (GTDB) release on
Uppmax, do:

.. code-block:: bash

    mkdir -p resources/gtdb
    ln -s /sw/data/GTDB/R202/$CLUSTER/release202/* resources/gtdb/


Then make sure your config file contains:

.. code-block:: bash

    binning:
        gtdbtk: True

nr
===================

Uppmax provides monthly snapshots of the ``nr`` non-redundant database. While
the formatted file cannot be used directly with the ``nbis-meta`` workflow you
can save time by making use of the already downloaded fasta file.

To use the latest snapshot of ``nr`` for taxonomic annotation of contigs, do:

.. code-block:: bash

    mkdir resources/nr
    ln -s /sw/data/diamond_databases/Blast/latest/download/nr.gz resources/nr/nr.fasta.gz
    file /sw/data/diamond_databases/Blast/latest | egrep -o "[0-9]{8}-[0-9]{6}" > resources/nr/nr.version

Then update the ``taxonomy`` section in your config file to use the ``nr`` database:

.. code-block:: yaml

    taxonomy:
        database: "nr"


Uniref90
========

The UniRef90 database is clustered at 90% sequence identity and Uppmax provides
downloaded fasta files that can be used directly with the workflow.

To use the latest snapshot of ``UniRef90`` for taxonomic annotation of contigs,
do:

.. code-block:: bash

    mkdir resources/uniref90
    ln -s /sw/data/diamond_databases/UniRef90/latest/download/uniref90.fasta.gz resources/uniref90/uniref90.fasta.gz
    file /sw/data/diamond_databases/UniRef90/latest | egrep -o "[0-9]{4}_[0-9]{2}" > resources/uniref90/uniref90.version

Then update the ``taxonomy`` section in your config file to use the ``uniref90``
database:

.. code-block:: yaml

    taxonomy:
        database: "uniref90"
