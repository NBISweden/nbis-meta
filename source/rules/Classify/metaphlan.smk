localrules:
    merge_metaphlan

rule metaphlan_pe:
    input:
        R1 = opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R1"+PREPROCESS+".fastq.gz"),
        R2 = opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R2"+PREPROCESS+".fastq.gz"),
        db = expand(opj(config["resource_path"], "metaphlan", "{index}.{s}.bt2"),
                    index = config["metaphlan_index"],
                    s = ["1","2","3","4","rev.1","rev.2"])
    output:
        tsv = opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_pe.tsv"),
        bt2 = opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_pe.bt2")
    log:
        opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_pe.log")
    params:
        dir = opj(config["resource_path"], "metaphlan")
    conda:
        "../../../envs/metaphlan.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metaphlan {input.R1},{input.R2} --bowtie2db {params.dir} \
            --nproc {threads} --input_type fastq -o {output.tsv} \
             --bowtie2out {output.bt2} > {log} 2>&1
        """

rule metaphlan_se:
    input:
        se = opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_se"+PREPROCESS+".fastq.gz"),
        db = expand(opj(config["resource_path"], "metaphlan", "{index}.{s}.bt2"),
                    index = config["metaphlan_index"],
                    s = ["1","2","3","4","rev.1","rev.2"])
    output:
        tsv = opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_se.tsv")
    log:
        opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_se.log")
    params:
        dir = opj(config["resource_path"], "metaphlan")
    conda:
        "../../../envs/metaphlan.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """       
        metaphlan {input.se} --bowtie2db {params.dir} \
            --nproc {threads} --input_type fastq -o {output.tsv} > {log} 2>&1
        """

rule merge_metaphlan:
    input:
        get_all_files(samples, opj(config["results_path"], "metaphlan"),
                      suffix=".tsv", nested=True)
    output:
        opj(config["report_path"], "metaphlan", "metaphlan.tsv")
    conda:
        "../../../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """