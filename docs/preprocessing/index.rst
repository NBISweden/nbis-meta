Preprocessing reads
===================

Read preprocessing can be run as a single part of the workflow using the command::

    snakemake --configfile <yourconfigfile> preprocess


Input reads can be trimmed using either `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
or `Cutadapt <https://github.com/marcelm/cutadapt>`_.
 
Trimmomatic
-----------
The settings specific to Trimmomatic are:

:code:`trimmomatic:` Set to :code:`True` to preprocess reads using Trimmomatic.

:code:`trimmomatic_home:` The directory where Trimmomatic stores the .jar file and adapter sequences. If you don't know it
leave it blank to let the pipeline attempt to locate it.

:code:`trim_adapters:` Set to 'True' to perform adapter trimming.

:code:`pe_adapter_params:` The adapter trim settings for paired end reads. This is what follows the 'ILLUMINACLIP' flag.
The default "2:30:15" will look for seeds with a maximum of **2** mismatches, and clip if extended seeds reach a score
of **30** for paired-end reads or **10**.

:code:`pe_pre_adapter_params:` Trim settings to be performed prior to adapter trimming. See the
`Trimmomatic manual <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf>`_ for
possible settings. As an example, to trim the first 10 bp from the start of reads set this to :code:`HEADCROP:10`.

:code:`pe_post_adapter_params:` Trim settings to be performed after adapter trimming. To for instance set a 50 bp
threshold on the minimum lenghts of reads after all trimming is done, set this to :code:`MINLEN:50`.

The :code:`se_adapter_params:` and :code:`se_post_adapter_params:` settings are the same as above but for single-end reads.

Cutadapt
--------

:code:`cutadapt:` Set to :code:`True` to run preprocessing with cutadapt.

.. Note:: Trimmomatic has priority in the preprocessing so if both Trimmomatic and Cutadapt are set to True, only Trimmomatic will be run.

:code:`adapter_sequence:` Adapter sequence for trimming. By default the workflow uses the Illumina TruSeq Universal Adapter.

:code:`rev_adapter_sequence:` 3' adapter to be removed from second read in a pair.

:code:`cutadapt_error_rate:` Maximum allowed error rate as value between 0 and 1. Defaults to 0.1. Increasing this value removes more adapters.

Phix filtering
--------------
:code:`phix_filter`: Set to True to filter out sequences mapping to the PhiX genome.

Fastuniq
--------
:code:`fastuniq:` Set to True to run de-duplication of paired reads using Fastuniq.

.. Note:: Fastuniq only runs with paired-end reads so if your data contains single-end samples the sequences will just be propagated downstream without Fastuniq processing.

SortMeRNA
---------
`SortMeRNA <https://github.com/biocore/sortmerna>`_ finds rRNA reads by
aligning to several rRNA databases. It can output aligning, rRNA, reads
and non-aligning, non_rRNA, reads to different output files allowing
you to filter your sequences.

:code:`sortmerna:` Set to :code:`True` to filter your raw reads with SortMeRNA.

:code:`sortmerna_keep:` Sortmerna produces files with reads aligning to rRNA ('rRNA' extension)
and not aligning to rRNA ('non_rRNA') extension. With the :code:`sortmerna_keep` setting
you specify which set of sequences you want to use for downstream analyses
('non_rRNA' or 'rRNA')

:code:`sortmerna_remove_filtered:` Set to True to remove the filtered reads (i.e. the reads NOT specified
in 'keep:')

:code:`sortmerna_dbs:` Databases to use for rRNA filtering. Can include:

- rfam-5s-database-id98.fasta
- rfam-5.8s-database-id98.fasta
- silva-arc-16s-id95.fasta
- silva-arc-23s-id98.fasta
- silva-bac-16s-id90.fasta
- silva-bac-23s-id98.fasta
- silva-euk-18s-id95.fasta
- silva-euk-28s-id98.fasta

:code:`sortmerna_paired_strategy:` How to handle read-pairs where mates are classified differently. If set to
:code:`paired_in` both reads in a pair are put into the 'rRNA' bin if one of them aligns (i.e. more strict)
while :code:`paired_out` puts both reads in the 'other' bin.

:code:`sortmerna_params:` Extra parameters to use for the sortmerna step.

Markduplicates
--------------
This is technically post-processing but to remove duplicates prior to producing read counts of ORFs called on assembled
contigs you can set :code:`markduplicates:True`. The :code:`picard_jar` and :code:`picard_path` can most often be left blank as the workflow will automatically identify
these paths in your conda environment. However, if you run into trouble with this step try searching for
:code:`picard.jar` and its directory.