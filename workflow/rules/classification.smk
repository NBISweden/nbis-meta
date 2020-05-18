from scripts.common import classify_input, get_kraken_index_url, get_krona_input

localrules:
    download_kraken_build,
    classifier2krona,
    all2krona

##### classification master rule #####

rule classify:
    input:
        classify_input(config)

##### kraken #####

rule download_kraken_build:
    """Downloads pre-built kraken2 index"""
    output:
        expand(opj(config["resource_path"], "kraken", "prebuilt",
                   config["kraken_prebuilt"], "{n}.k2d"),
               n=["hash", "opts", "taxo"])
    log:
        opj(config["resource_path"], "kraken", "prebuilt",
                   config["kraken_prebuilt"], "download.log")
    params:
        dir=lambda w, output: os.path.dirname(output[0]),
        tar=opj(config["scratch_path"],
                "{base}.tgz".format(base=config["kraken_prebuilt"])),
        url=get_kraken_index_url(config["kraken_prebuilt"]),
        db_version=get_kraken_index_url(config["kraken_prebuilt"], version=True),
        tmpdir = opj(config["scratch_path"], "kraken_db")
    shell:
         """
         mkdir -p {params.tmpdir}
         curl -L -o {params.tar} {params.url} > {log} 2>&1
         tar -C {params.tmpdir} -xf {params.tar}
         mv {params.tmpdir}/*/* {params.dir}
         rm -r {params.tar} {params.tmpdir}
         echo {params.db_version} > {params.dir}/version
         """

rule kraken_build_standard:
    output:
        expand(opj(config["resource_path"], "kraken", "standard", "{n}.k2d"),
               n=["hash", "opts", "taxo"])
    log:
        build=opj(config["resource_path"], "kraken", "standard", "build.log"),
        clean=opj(config["resource_path"], "kraken", "standard", "clean.log"),
    params:
        dir=lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/kraken.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*24
    shell:
        """
        kraken2-build --standard --db {params.dir} --threads {threads} > {log.build} 2>&1
        kraken2-build --clean {params.dir} > {log.clean} 2>&1
        """

rule kraken_pe:
    input:
        R1=opj(config["intermediate_path"], "preprocess",
               "{sample}_{run}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["intermediate_path"], "preprocess",
               "{sample}_{run}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["kraken_index_path"], "{n}.k2d"),
                  n=["hash", "opts", "taxo"])
    output:
        opj(config["results_path"], "kraken", "{sample}_{run}_pe.out"),
        opj(config["results_path"], "kraken", "{sample}_{run}_pe.kreport")
    log:
        opj(config["results_path"], "kraken", "{sample}_{run}_pe.log")
    params:
        db=opj(config["kraken_index_path"]),
        mem=config["kraken_params"]
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/kraken.yml"
    shell:
        """
        kraken2 {params.mem} --db {params.db} --output {output[0]} \
            --report {output[1]} --gzip-compressed \
            --threads {threads} --paired {input.R1} {input.R2} > {log} 2>&1
        """

rule kraken_se:
    input:
        se=opj(config["intermediate_path"], "preprocess",
               "{sample}_{run}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["kraken_index_path"], "{n}.k2d"),
                  n=["hash", "opts", "taxo"])
    output:
        opj(config["results_path"], "kraken", "{sample}_{run}_se.out"),
        opj(config["results_path"], "kraken", "{sample}_{run}_se.kreport")
    log:
        opj(config["results_path"], "kraken", "{sample}_{run}_se.log")
    params:
        db=opj(config["kraken_index_path"]),
        mem=config["kraken_params"]
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/kraken.yml"
    shell:
        """
        kraken2 {params.mem} --db {params.db} --output {output[0]} \
            --report {output[1]} --gzip-compressed \
            --threads {threads} {input.se} > {log} 2>&1
        """


##### krona #####

rule krona_taxonomy:
    output:
        tab = opj(config["resource_path"], "krona", "taxonomy.tab")
    log:
        opj(config["resource_path"], "krona", "taxonomy.log")
    params:
        taxdir = lambda w, output: os.path.dirname(output.tab)
    conda:
        "../envs/krona.yml"
    singularity:
        "docker://continuumio/miniconda3:4.8.2"
    shell:
        """
        ktUpdateTaxonomy.sh {params.taxdir} >{log} 2>&1
        """

rule classifier2krona:
    input:
        opj(config["results_path"], "{classifier}",
            "{sample}_{run}_{seq_type}.kreport"),
        opj("resources", "krona", "taxonomy.tab")
    output:
        opj(config["results_path"], "{classifier}",
            "{sample}_{run}_{seq_type}.html")
    params:
        tax="resources/krona"
    conda:
        "../envs/krona.yml"
    singularity:
        "docker://continuumio/miniconda3:4.8.2"
    shell:
        """
        ktImportTaxonomy -t 5 -m 3 \
            -tax {params.tax} -o {output[0]} \
            {input[0]},{wildcards.sample}_{wildcards.run}
        """

rule all2krona:
    input:
        f=get_all_files(samples, opj(config["results_path"],
                                    "{classifier}"), ".kreport"),
        h=get_all_files(samples, opj(config["results_path"],
                                    "{classifier}"), ".html"),
        t=opj("resources", "krona", "taxonomy.tab")
    output:
        opj(config["report_path"], "{classifier}", "{classifier}.krona.html")
    log:
        opj(config["report_path"], "{classifier}", "{classifier}.krona.log")
    params:
        tax="resources/krona",
        input_string=get_krona_input(config, samples, "{classifier}")
    conda:
        "../envs/krona.yml"
    singularity:
        "docker://continuumio/miniconda3:4.8.2"
    shell:
         """
         ktImportTaxonomy \
            -t 5 -m 3 -tax {params.tax} -o {output[0]} \
            {params.input_string} > {log} 2>&1
         """
