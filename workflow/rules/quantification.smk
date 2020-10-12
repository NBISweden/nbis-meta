localrules:
    quantify,
    write_featurefile,
    aggregate_featurecount,
    clean_featurecount,
    count_features,
    normalize_features,
    sum_to_taxa,
    sum_to_rgi

##### quantify master rule #####

rule quantify:
    input:
        expand(opj(config["paths"]["results"], "annotation", "{assembly}",
                   "gene_{counts_type}.tsv"),
               assembly=assemblies.keys(), counts_type=["counts", "rpkm"])


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

rule featurecount:
    input:
        gff=opj(config["paths"]["results"], "annotation", "{assembly}",
                "final_contigs.features.gff"),
        bam=opj(config["paths"]["results"], "assembly", "{assembly}",
                "mapping", "{sample}_{unit}_{seq_type}"+POSTPROCESS+".bam")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.tsv"),
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.tsv.summary")
    log:
        opj(config["paths"]["results"], "assembly", "{assembly}",
            "mapping", "{sample}_{unit}_{seq_type}.fc.log")
    threads: 4
    params:
        tmpdir=config["paths"]["temp"],
        setting=lambda wildcards: "-B -p" if wildcards.seq_type == "pe" else ""
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
    conda:
        "../envs/quantify.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        featureCounts -a {input.gff} -o {output[0]} -t CDS -g gene_id -M \
            {params.setting} -T {threads} --tmpDir {params.tmpdir} {input.bam} > {log} 2>&1
        """

rule clean_featurecount:
    input:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.tsv")
    output:
        opj(config["paths"]["results"], "assembly", "{assembly}", "mapping",
            "{sample}_{unit}_{seq_type}.fc.clean.tsv")
    script:
        "../scripts/quantification_utils.py"

rule aggregate_featurecount:
    """Aggregates all cleaned count files from featureCounts"""
    input:
        get_all_files(samples=samples,
                      dir=opj(config["paths"]["results"],"assembly", "{assembly}", "mapping"),
                      suffix=".fc.clean.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "gene_counts.tsv")
    script:
        "../scripts/quantification_utils.py"

rule rpkm:
    """
    Calculate RPKM for genes in an assembly
    """
    input:
        opj(config["paths"]["results"], "annotation", "{assembly}", "gene_counts.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "gene_rpkm.tsv")
    log:
        opj(config["paths"]["results"], "annotation", "{assembly}", "rpkm.log")
    params:
        method = "RPKM"
    conda:
        "../envs/normalize.yml"
    script:
        "../scripts/normalize.R"

rule count_features:
    """
    Sums read counts for gene annotation features such as pfam, KOs etc.
    """
    input:
        abund=opj(config["paths"]["results"], "annotation", "{assembly}", "gene_counts.tsv"),
        annot=opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.counts.tsv")
    script:
        "../scripts/quantification_utils.py"

rule normalize_features:
    """
    Normalizes counts of features using TMM, REL or CSS
    """
    input:
        opj(config["paths"]["results"], "annotation", "{assembly}", "{db}.parsed.counts.tsv")
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
    """
    Sums read counts and RPKM values for genes to assigned taxonomy
    """
    input:
        tax=opj(config["paths"]["results"], "annotation", "{assembly}", "taxonomy",
            "orfs.{db}.taxonomy.tsv".format(db=config["taxonomy"]["database"])),
        abund=opj(config["paths"]["results"], "annotation", "{assembly}", "gene_{counts_type}.tsv")
    output:
        opj(config["paths"]["results"], "annotation", "{assembly}", "taxonomy", "tax.{counts_type}.tsv")
    script:
        "../scripts/quantification_utils.py"
