localrules:
    quantify,
    write_featurefile,
    samtools_stats,
    normalize_featurecount,
    aggregate_featurecount,
    sum_to_taxa,
    quantify_features,
    sum_to_rgi

##### quantify master rule #####

rule quantify:
    input:
        expand(opj(config["paths"]["results"], "annotation", "{assembly}",
                   "fc.{fc_type}.tsv"),
               assembly=assemblies.keys(), fc_type=["tpm", "raw"])


rule write_featurefile:
    input:
        opj(config["paths"]["results"], "annotation", "{assembly}",
            "final_contigs.gff")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}",
            "final_contigs.features.gff")
    script:
        "../scripts/quantification_utils.py"

##### markduplicates #####

rule remove_mark_duplicates:
    input:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping", "{sample}_{unit}_{seq_type}.bam")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping", "{sample}_{unit}_{seq_type}.markdup.bam"),
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping", "{sample}_{unit}_{seq_type}.markdup.bam.bai"),
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping", "{sample}_{unit}_{seq_type}.markdup.metrics")
    log:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping", "{sample}_{unit}_{seq_type}.markdup.log")
    params:
        header=opj(config["paths"]["temp"], "{assembly}", "{sample}_{unit}_{seq_type}.header"),
        rehead_bam=opj(config["paths"]["temp"], "{assembly}", "{sample}_{unit}_{seq_type}.rehead.bam"),
        temp_bam=opj(config["paths"]["temp"], "{assembly}", "{sample}_{unit}_{seq_type}.markdup.bam"),
        temp_sort_bam=opj(config["paths"]["temp"], "{assembly}", "{sample}_{unit}_{seq_type}.markdup.re_sort.bam"),
        temp_dir=opj(config["paths"]["temp"], "{assembly}")
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/quantify.yml"
    shell:
        """
        mkdir -p {params.temp_dir}
        # Fix bam header
        samtools view -H {input} | egrep -v "^@PG" > {params.header}
        samtools reheader -P {params.header} {input} > {params.rehead_bam}
        # Set memory max
        mem="-Xmx$((6 * {threads}))g"
        java -Xms2g $mem -XX:ParallelGCThreads={threads} \
            -jar $CONDA_PREFIX/share/picard-*/picard.jar MarkDuplicates \
            I={params.rehead_bam} M={output[2]} O={params.temp_bam} REMOVE_DUPLICATES=TRUE \
            USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE ASSUME_SORT_ORDER=coordinate \
            PROGRAM_RECORD_ID=null ADD_PG_TAG_TO_READS=FALSE 2> {log}
        # Re sort the bam file using samtools
        samtools_threads="$(({threads} - 1))"
        samtools sort -@ $samtools_threads -o {params.temp_sort_bam} {params.temp_bam} > /dev/null 2>&1
        # Index the bam file
        samtools index {params.temp_sort_bam}
        mv {params.temp_sort_bam} {output[0]}
        mv {params.temp_sort_bam}.bai {output[1]}
        rm {params.temp_bam} {params.rehead_bam} {params.header}
        """

##### featurecounts #####

rule featurecount_pe:
    input:
        gff=opj(config["paths"]["results"], "annotation", "{assembly}",
                "final_contigs.features.gff"),
        bam=opj(config["paths"]["results"], "assembly", "{assembly}",
                "mapping", "{sample}_{unit}_pe"+POSTPROCESS+".bam")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_pe.fc.tsv"),
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_pe.fc.tsv.summary")
    log:
        opj(config["paths"]["results"], "assembly", "{assembly}",
            "mapping", "{sample}_{unit}_pe.fc.log")
    threads: 4
    params: tmpdir=config["paths"]["temp"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
    conda:
        "../envs/quantify.yml"
    shell:
        """
        featureCounts -a {input.gff} -o {output[0]} -t CDS -g gene_id -M -p \
            -B -T {threads} --tmpDir {params.tmpdir} {input.bam} > {log} 2>&1
        """

rule featurecount_se:
    input:
        gff=opj(config["paths"]["results"], "annotation", "{assembly}",
                "final_contigs.features.gff"),
        bam=opj(config["paths"]["results"], "assembly", "{assembly}",
                "mapping", "{sample}_{unit}_se"+POSTPROCESS+".bam")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}",
            "mapping", "{sample}_{unit}_se.fc.tsv"),
        opj(config["paths"]["results"], "assembly", "{assembly}",
            "mapping", "{sample}_{unit}_se.fc.tsv.summary")
    log:
        opj(config["paths"]["results"], "assembly", "{assembly}",
            "mapping", "{sample}_{unit}_se.fc.log")
    threads: 4
    params: tmpdir=config["paths"]["temp"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
    conda:
        "../envs/quantify.yml"
    shell:
        """
        featureCounts -a {input.gff} -o {output[0]} -t CDS -g gene_id \
            -M -T {threads} --tmpDir {params.tmpdir} {input.bam} > {log} 2>&1
        """

rule samtools_stats:
    input:
        opj(config["paths"]["results"], "assembly", "{assembly}",
                "mapping", "{sample}_{unit}_{seq_type}"+POSTPROCESS+".bam")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}",
                "mapping", "{sample}_{unit}_{seq_type}"+POSTPROCESS+".bam.stats")
    conda:
        "../envs/quantify.yml"
    shell:
        """
        samtools stats {input} > {output}
        """

rule normalize_featurecount:
    input:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.tsv"),
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}"+POSTPROCESS+".bam.stats")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.tpm.tsv"),
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.raw.tsv")
    log:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.norm.log")
    script:
        "../scripts/quantification_utils.py"

rule aggregate_featurecount:
    """Aggregates raw and normalized featureCounts files"""
    input:
        raw_files=get_all_files(samples, opj(config["paths"]["results"],
                                             "assembly", "{assembly}", "mapping"),
                                ".fc.raw.tsv"),
        tpm_files=get_all_files(samples, opj(config["paths"]["results"],
                                             "assembly", "{assembly}", "mapping"),
                                ".fc.tpm.tsv"),
        gff_file=opj(config["paths"]["results"], "annotation", "{assembly}",
                     "final_contigs.features.gff")
    output:
        raw=opj(config["paths"]["results"], "annotation", "{assembly}", "fc.raw.tsv"),
        tpm=opj(config["paths"]["results"], "annotation", "{assembly}", "fc.tpm.tsv")
    script:
        "../scripts/quantification_utils.py"

rule quantify_features:
    input:
        abund=opj(config["paths"]["results"], "annotation", "{assembly}", "fc.{fc_type}.tsv"),
        annot=opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.{fc_type}.tsv")
    shell:
        """
        python workflow/scripts/eggnog-parser.py \
            quantify {input.abund} {input.annot} {output[0]}
        """

rule normalize:
    input:
        opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.raw.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.{norm_method}.tsv")
    log:
        opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.{norm_method}.log")
    params:
        method = "{norm_method}"
    conda:
        "../envs/normalize.yml"
    script:
        "../scripts/normalize.R"

rule sum_to_taxa:
    input:
        tax=opj(config["paths"]["results"], "annotation", "{assembly}", "taxonomy",
            "orfs.{db}.taxonomy.tsv".format(db=config["taxonomy"]["database"])),
        abund=opj(config["paths"]["results"], "annotation", "{assembly}", "fc.{fc_type}.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "taxonomy", "tax.{fc_type}.tsv")
    script:
        "../scripts/quantification_utils.py"

rule sum_to_rgi:
    input:
        annot=opj(config["paths"]["results"], "annotation", "{assembly}", "rgi.out.txt"),
        abund=opj(config["paths"]["results"], "annotation", "{assembly}", "fc.{fc_type}.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "rgi.{fc_type}.tsv")
    script:
        "../scripts/quantification_utils.py"

