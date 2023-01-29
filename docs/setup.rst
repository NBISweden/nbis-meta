###############
Getting started
###############

.. _config-file:

**********************
The configuration file
**********************

The workflow comes with a configuration file ``config/config.yaml``. You can
either modify this file directly or create a new config file and point snakemake
to it with ``snakemake --configfile <your-config-file>.yaml``. The configuration
is `validated <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation>`_
with a schema found under ``workflow/schemas/config.schema.yaml`` which is also
used to set default values. Specifying a configfile with ``--configfile`` extends
and overwrites the default settings so that any configfile you specify only
needs to contain the parameters you want to change.

For example, the default config file begins with:

.. code-block:: yaml

    sample_list: config/samples.tsv

    paths:
      results: "results"

A new config file, named *e.g.* ``conf_test.yaml``, and containing **only** the
lines:

.. code-block:: yaml

    paths:
        results: "test"

will then use ``test/`` as output directory for results when the workflow is run
as:

.. code-block:: bash

    snakemake --use-conda --configfile conf_test.yaml -j 4

.. _sample-list:

***************
The sample list
***************

The workflow requires you to supply a file with some information on your
samples. The path to this file is specified with the ``sample_list:`` parameter
in the configuration file.

As a minimum the file must contain the columns ``sample``, ``unit`` and ``fq1``.
Paired-end samples also require the ``fq2`` column for the mate 2 read file.

A very basic sample file may look like this:

+---------+------+----+--------------------------------+-------------------------------------+
| sample  | unit | fq1                                 |    fq2                              |
+=========+======+=====================================+=====================================+
| sample1 |  1   | examples/data/sample1_1_R1.fastq.gz | examples/data/sample1_1_R2.fastq.gz |
+---------+------+-------------------------------------+-------------------------------------+
| sample2 |  1   | examples/data/sample2_11_R1.fastq.gz| examples/data/sample2_11_R2.fastq.gz|
+---------+------+-------------------------------------+-------------------------------------+
| sample3 |  1   | examples/data/sample3_21_R1.fastq.gz|                                     |
+---------+------+-------------------------------------+-------------------------------------+

The columns ``fq1`` and ``fq2`` (for paired-end data) specify the paths to fastq
files which will be used as input to the workflow.

Specifying assemblies
=====================

You may also include a column named ``assembly`` with names for assemblies to
create. The assembly field can be comma-separated entries specifying several
assemblies per sample/unit combination.

For instance, with a sample file like this:

+---------+------+-------------+--------------------------------------+-------------------------------------+
| sample  | unit | assembly    | fq1                                  |    fq2                              |
+=========+======+=============+======================================+=====================================+
| sample1 |  1   | sample1,all | examples/data/sample1_1_R1.fastq.gz  | examples/data/sample1_1_R2.fastq.gz |
+---------+------+-------------+--------------------------------------+-------------------------------------+
| sample2 |  1   | sample2,all | examples/data/sample2_11_R1.fastq.gz | examples/data/sample2_11_R2.fastq.gz|
+---------+------+-------------+--------------------------------------+-------------------------------------+
| sample3 |  1   | sample3,all | examples/data/sample3_21_R1.fastq.gz |                                     |
+---------+------+-------------+--------------------------------------+-------------------------------------+

a total of 4 assemblies will be generated:
- ``sample1``, ``sample2`` and ``sample3`` with input only from each sample respectively, and
- ``all`` with input from all samples
