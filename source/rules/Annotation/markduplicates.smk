rule remove_mark_duplicates:
    input:
        opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_{seq_type}.bam")
    output:
        opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_{seq_type}.markdup.bam"),
        opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_{seq_type}.markdup.bam.bai"),
        opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_{seq_type}.markdup.metrics")
    log:
        opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_{seq_type}.markdup.log")
    params:
        temp_bam=opj(config["tmpdir"],"{group}","{sample}_{run}_{seq_type}.markdup.bam"),
        temp_sort_bam=opj(config["tmpdir"],"{group}", "{sample}_{run}_{seq_type}.markdup.re_sort.bam"),
        temp_dir=opj(config["tmpdir"],"{group}")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/quantify.yml"
    threads: 10
    shell:
        """
        mkdir -p {params.temp_dir}
        java \
            -Xms2g -Xmx64g \
            -XX:ParallelGCThreads={threads} \
            -jar $CONDA_PREFIX/share/picard-*/picard.jar \
            MarkDuplicates \
            I={input} \
            M={output[2]} \
            O={params.temp_bam} \
            REMOVE_DUPLICATES=TRUE \
            USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
            ASO=coordinate \
            2> {log}
        # Re sort the bam file using samtools
        samtools sort -o {params.temp_sort_bam} {params.temp_bam}
        # Index the bam file
        samtools index {params.temp_sort_bam}
        mv {params.temp_sort_bam} {output[0]}
        mv {params.temp_sort_bam}.bai {output[1]}
        rm {params.temp_bam}
        """