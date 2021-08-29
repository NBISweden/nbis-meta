localrules:
    quantify,
    write_featurefile,
    clean_featurecount,
    aggregate_featurecount,
    rpkm,
    count_features,
    edger_normalize_features,
    css_normalize_features,
    sum_to_taxa

##### quantify master rule #####

rule quantify:
    input:
        expand("{results}/annotation/{assembly}/gene_{counts_type}.tsv",
            results=[config["paths"]["results"]], assembly=assemblies.keys(),
            counts_type=["counts", "rpkm"])


rule write_featurefile:
    input:
        results+"/annotation/{assembly}/final_contigs.gff"
    output:
        results+"/annotation/{assembly}/final_contigs.features.gff"
    script:
        "../scripts/quantification_utils.py"

##### markduplicates #####

rule remove_mark_duplicates:
    input:
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.bam"
    output:
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.markdup.bam",
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.markdup.bam.bai",
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.markdup.metrics"
    log:
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.markdup.log"
    params:
        header=temppath+"/{assembly}/{sample}_{unit}_{seq_type}.header",
        rehead_bam=temppath+"/{assembly}/{sample}_{unit}_{seq_type}.rehead.bam",
        temp_bam=temppath+"/{assembly}/{sample}_{unit}_{seq_type}.markdup.bam",
        temp_sort_bam=temppath+"/{assembly}/{sample}_{unit}_{seq_type}.markdup.re_sort.bam",
        temp_dir=temppath+"/{assembly}"
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
        gff=results+"/annotation/{assembly}/final_contigs.features.gff",
        bam=results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}"+POSTPROCESS+".bam"
    output:
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.fc.tsv",
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.fc.tsv.summary"
    log:
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.fc.log"
    threads: 4
    params:
        tmpdir=config["paths"]["temp"],
        setting=lambda wildcards: "-Q 10 -B -p" if wildcards.seq_type == "pe" else ""
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
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.fc.tsv"
    output:
        results+"/assembly/{assembly}/mapping/{sample}_{unit}_{seq_type}.fc.clean.tsv"
    script:
        "../scripts/quantification_utils.py"

rule aggregate_featurecount:
    """Aggregates all cleaned count files from featureCounts"""
    input:
        get_all_files(samples=samples,
                      directory=results+"/assembly/{assembly}/mapping",
                      suffix=".fc.clean.tsv")
    output:
        results+"/annotation/{assembly}/gene_counts.tsv"
    script:
        "../scripts/quantification_utils.py"

rule rpkm:
    """
    Calculate RPKM for genes in an assembly
    """
    input:
        results+"/annotation/{assembly}/gene_counts.tsv"
    output:
        results+"/annotation/{assembly}/gene_rpkm.tsv"
    log:
        results+"/annotation/{assembly}/rpkm.log"
    params:
        method = "RPKM"
    conda:
        "../envs/edger.yml"
    script:
        "../scripts/edger.R"

rule count_features:
    """
    Sums read counts for gene annotation features such as pfam, KOs etc.
    """
    input:
        abund=results+"/annotation/{assembly}/gene_counts.tsv",
        annot=results+"/annotation/{assembly}/{db}.parsed.tsv"
    output:
        results+"/annotation/{assembly}/{db}.parsed.counts.tsv"
    script:
        "../scripts/quantification_utils.py"

rule edger_normalize_features:
    """
    Normalizes counts of features using TMM and REL
    """
    input:
        results+"/annotation/{assembly}/{db}.parsed.counts.tsv"
    output:
        results+"/annotation/{assembly}/{db}.parsed.{norm_method}.tsv"
    log:
        results+"/annotation/{assembly}/{db}.parsed.{norm_method}.log"
    params:
        method = "{norm_method}"
    conda:
        "../envs/edger.yml"
    script:
        "../scripts/edger.R"

rule css_normalize_features:
    """
    Normalizes counts of features using CSS from metagenomeSeq
    """
    input:
        results+"/annotation/{assembly}/{db}.parsed.counts.tsv"
    output:
        results+"/annotation/{assembly}/{db}.parsed.CSS.tsv"
    log:
        results+"/annotation/{assembly}/{db}.parsed.CSS.log"
    conda:
        "../envs/metagenomeseq.yml"
    script:
        "../scripts/metagenomeseq.R"

rule sum_to_taxa:
    """
    Sums read counts and RPKM values for genes to assigned taxonomy
    """
    input:
        tax=expand(results+"/annotation/{{assembly}}/taxonomy/orfs.{db}.taxonomy.tsv",
            db=config["taxonomy"]["database"]),
        abund=results+"/annotation/{assembly}/gene_{counts_type}.tsv"
    output:
        results+"/annotation/{assembly}/taxonomy/tax.{counts_type}.tsv"
    script:
        "../scripts/quantification_utils.py"
