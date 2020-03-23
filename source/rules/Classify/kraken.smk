################
## Run Kraken ##
################
rule kraken_pe:
    input:
        R1=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["kraken_index_path"],"{n}.k2d"),
                  n=["hash","opts","taxo"])
    output:
        opj(config["results_path"],"kraken","{sample}_{run}_pe.out"),
        opj(config["results_path"],"kraken","{sample}_{run}_pe.kreport")
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    params:
        db=opj(config["kraken_index_path"]),
        mem=config["kraken_params"]
    conda:
        "../../../envs/kraken.yml"
    shell:
        """
        kraken2 \
            {params.mem} --db {params.db} --output {output[0]} \
            --report {output[1]} --gzip-compressed \
            --threads {threads} --paired {input.R1} {input.R2}
        """

rule kraken_se:
    input:
        se=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["kraken_index_path"],"{n}.k2d"),
                  n=["hash","opts","taxo"])
    output:
        opj(config["results_path"],"kraken","{sample}_{run}_se.out"),
        opj(config["results_path"],"kraken","{sample}_{run}_se.kreport")
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    params:
        db=opj(config["kraken_index_path"]),
        mem=config["kraken_params"]
    conda:
        "../../../envs/kraken.yml"
    shell:
        """
        kraken2 \
            {params.mem} --db {params.db} --output {output[0]} \
            --report {output[1]} --gzip-compressed \
            --threads {threads} {input.se}
        """