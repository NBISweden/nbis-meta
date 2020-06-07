from scripts.common import classify_input, krona_input, metaphlan_krona_string
from scripts.common import get_kraken_index_url, get_centrifuge_index_url

localrules:
    classify,
    download_kraken_build,
    download_centrifuge_build,
    centrifuge_kreport,
    merge_metaphlan,
    metaphlan2krona_table,
    metaphlan2krona,
    plot_metaphlan,
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
        expand(opj("resources", "kraken", "prebuilt",
                   config["kraken"]["prebuilt"], "{n}.k2d"),
               n=["hash", "opts", "taxo"])
    log:
        opj("resources", "kraken", "prebuilt",
                   config["kraken"]["prebuilt"], "download.log")
    params:
        dir=lambda w, output: os.path.dirname(output[0]),
        tar=opj(config["paths"]["temp"],
                "{base}.tgz".format(base=config["kraken"]["prebuilt"])),
        url=get_kraken_index_url(config["kraken"]["prebuilt"]),
        db_version=get_kraken_index_url(config["kraken"]["prebuilt"], version=True),
        tmpdir = opj(config["paths"]["temp"], "kraken_db")
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
        expand(opj("resources", "kraken", "standard", "{n}.k2d"),
               n=["hash", "opts", "taxo"])
    log:
        build=opj("resources", "kraken", "standard", "build.log"),
        clean=opj("resources", "kraken", "standard", "clean.log"),
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
        R1=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["kraken"]["index_path"], "{n}.k2d"),
                  n=["hash", "opts", "taxo"])
    output:
        opj(config["paths"]["results"], "kraken", "{sample}_{unit}_pe.out"),
        opj(config["paths"]["results"], "kraken", "{sample}_{unit}_pe.kreport")
    log:
        opj(config["paths"]["results"], "kraken", "{sample}_{unit}_pe.log")
    params:
        db=opj(config["kraken"]["index_path"]),
        mem=config["kraken"]["mem"]
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
        se=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["kraken"]["index_path"], "{n}.k2d"),
                  n=["hash", "opts", "taxo"])
    output:
        opj(config["paths"]["results"], "kraken", "{sample}_{unit}_se.out"),
        opj(config["paths"]["results"], "kraken", "{sample}_{unit}_se.kreport")
    log:
        opj(config["paths"]["results"], "kraken", "{sample}_{unit}_se.log")
    params:
        db=opj(config["kraken"]["index_path"]),
        mem=config["kraken"]["mem"]
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

##### centrifuge #####

rule download_centrifuge_build:
    """Downloads pre-built centrifuge index"""
    output:
        db=expand(opj(config["centrifuge"]["dir"],
                      "{base}.{i}.cf"), i=[1, 2, 3],
                  base=config["centrifuge"]["base"])
    log:
        opj(config["centrifuge"]["dir"], "download.log")
    params:
        dir=config["centrifuge"]["dir"],
        tar=opj(config["centrifuge"]["dir"],
                "{base}.tar.gz".format(base=config["centrifuge"]["base"])),
        url=get_centrifuge_index_url(config)
    shell:
        """
        curl -o {params.tar} {params.url} > {log} 2>&1
        tar -C {params.dir} -xf {params.tar} >>{log} 2>&1
        rm {params.tar}
        """
    
rule centrifuge_pe:
    input:
        R1=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["centrifuge"]["dir"], "{base}.{i}.cf"),
                  i=[1, 2, 3], base=config["centrifuge"]["base"])
    output:
        opj(config["paths"]["results"], "centrifuge", "{sample}_{unit}_pe.out"),
        opj(config["paths"]["results"], "centrifuge", "{sample}_{unit}_pe.report")
    log:
        opj(config["paths"]["results"], "centrifuge", "{sample}_{unit}_pe.log")
    params:
        prefix=opj(config["centrifuge"]["dir"],
                   "{base}".format(base=config["centrifuge"]["base"])),
        tmp_out=opj(config["paths"]["temp"], "{sample}_{unit}_pe.out"),
        tmpdir=config["paths"]["temp"],
        tmp_report=opj(config["paths"]["temp"], "{sample}_{unit}_pe.report"),
        k=config["centrifuge"]["max_assignments"]
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    conda:
        "../envs/centrifuge.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        centrifuge -k {params.k} -x {params.prefix} -1 {input.R1} -2 {input.R2} \
            -S {params.tmp_out} -p {threads} --report-file {params.tmp_report} \
            > {log} 2>&1
        mv {params.tmp_out} {output[0]}
        mv {params.tmp_report} {output[1]}
        """

rule centrifuge_se:
    input:
        se=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["centrifuge"]["dir"], "{base}.{i}.cf"),
                  i=[1, 2, 3], base=config["centrifuge"]["base"])
    output:
        opj(config["paths"]["results"], "centrifuge", "{sample}_{unit}_se.out"),
        opj(config["paths"]["results"], "centrifuge", "{sample}_{unit}_se.report")
    log:
        opj(config["paths"]["results"], "centrifuge", "{sample}_{unit}_se.log")
    params:
        prefix=opj(config["centrifuge"]["dir"],
                   "{base}".format(base=config["centrifuge"]["base"])),
        tmp_out=opj(config["paths"]["temp"], "{sample}_{unit}_se.out"),
        tmpdir=config["paths"]["temp"],
        tmp_report=opj(config["paths"]["temp"], "{sample}_{unit}_se.report"),
        k=config["centrifuge"]["max_assignments"]
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    conda:
        "../envs/centrifuge.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        centrifuge -k {params.k} -U {input.se} -x {params.prefix} \
            -S {params.tmp_out} -p {threads} --report-file {params.tmp_report} \
            > {log} 2>&1
        mv {params.tmp_out} {output[0]}
        mv {params.tmp_report} {output[1]}
        """

rule centrifuge_kreport:
    input:
        f=opj(config["paths"]["results"], "centrifuge", 
              "{sample}_{unit}_{seq_type}.out"),
        db=expand(opj(config["centrifuge"]["dir"],
                      "{base}.{i}.cf"), 
                  i=[1, 2, 3], base=config["centrifuge"]["base"])
    output:
        opj(config["paths"]["results"], "centrifuge", 
            "{sample}_{unit}_{seq_type}.kreport")
    params:
        min_score=config["centrifuge"]["min_score"],
        prefix=opj(config["centrifuge"]["dir"],
                   "{base}".format(base=config["centrifuge"]["base"]))
    conda:
        "../envs/centrifuge.yml"
    shell:
        """
        centrifuge-kreport --min-score {params.min_score} -x {params.prefix} \
            {input.f} > {output[0]}
        """

##### metaphlan #####

rule build_metaphlan:
    """
    Download and build the metaphlan bowtie2 database
    """
    output:
        expand(opj("resources", "metaphlan", "{index}.{s}.bt2"),
               index=config["metaphlan"]["index"],
               s=["1", "2", "3", "4", "rev.1", "rev.2"])
    log:
        opj("resources", "metaphlan", "mpa.log")
    params:
        dir=opj("resources", "metaphlan"),
        index=config["metaphlan"]["index"]
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        metaphlan --install --bowtie2db {params.dir} \
            --nproc {threads} -x {params.index} >{log} 2>&1
        """

rule metaphlan_pe:
    input:
        R1=opj(config["paths"]["results"], "intermediate", "preprocess",
                 "{sample}_{unit}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["paths"]["results"], "intermediate", "preprocess",
                 "{sample}_{unit}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj("resources", "metaphlan", "{index}.{s}.bt2"),
                  index=config["metaphlan"]["index"],
                  s=["1", "2", "3", "4", "rev.1", "rev.2"])
    output:
        tsv=opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_pe.tsv"),
        bt2=opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_pe.bt2")
    log:
        opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_pe.log")
    params:
        dir=opj("resources", "metaphlan")
    conda:
        "../envs/metaphlan.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metaphlan {input.R1},{input.R2} --bowtie2db {params.dir} --add_viruses \
            --force --nproc {threads} --input_type fastq -o {output.tsv} \
             --bowtie2out {output.bt2} > {log} 2>&1
        """

rule metaphlan_se:
    input:
        se=opj(config["paths"]["results"], "intermediate", "preprocess",
               "{sample}_{unit}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj("resources", "metaphlan", "{index}.{s}.bt2"),
                  index = config["metaphlan"]["index"],
                  s = ["1", "2", "3", "4", "rev.1", "rev.2"])
    output:
        tsv=opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_se.tsv"),
        bt2=opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_se.bt2")
    log:
        opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_se.log")
    params:
        dir = opj("resources", "metaphlan")
    conda:
        "../envs/metaphlan.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """       
        metaphlan {input.se} --bowtie2db {params.dir} --add_viruses --force \
            --nproc {threads} --input_type fastq -o {output.tsv} \
             --bowtie2out {output.bt2} > {log} 2>&1
        """

rule merge_metaphlan:
    input:
        get_all_files(samples, opj(config["paths"]["results"], "metaphlan"), ".tsv")
    output:
        opj(config["paths"]["results"], "report", "metaphlan", "metaphlan.tsv")
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """

rule metaphlan2krona_table:
    input:
        opj(config["paths"]["results"], "metaphlan", "{sample}_{unit}_{seq_type}.tsv")
    output:
        temp(opj(config["paths"]["results"], "metaphlan",
                 "{sample}_{unit}_{seq_type}.krona"))
    script:
        "../scripts/classification_utils.py"

rule metaphlan2krona:
    input:
        files = get_all_files(samples, opj(config["paths"]["results"], "metaphlan"), ".krona"),
        db = opj("resources", "krona", "taxonomy.tab")
    output:
        opj(config["paths"]["results"], "report", "metaphlan", "metaphlan.html")
    log:
        opj(config["paths"]["results"], "report", "metaphlan", "krona.log")
    conda:
        "../envs/krona.yml"
    params:
        input_string = lambda w, input: metaphlan_krona_string(input.files),
        dbdir = lambda w, input: os.path.dirname(input.db)
    shell:
        """
        ktImportTaxonomy -t 1 -m 2 -o {output} -tax {params.dbdir} \
            {params.input_string} > {log} 2>&1
        """

rule plot_metaphlan:
    input:
        opj(config["paths"]["results"], "report", "metaphlan", "metaphlan.tsv")
    output:
        opj(config["paths"]["results"], "report", "metaphlan", "metaphlan.pdf")
    params:
        rank=config["metaphlan"]["plot_rank"]
    conda:
        "../envs/plotting.yml"
    notebook:
        "../notebooks/metaphlan.py.ipynb"
        
##### krona #####

rule krona_taxonomy:
    output:
        tab=opj("resources", "krona", "taxonomy.tab")
    log:
        opj("resources", "krona", "taxonomy.log")
    params:
        taxdir=lambda w, output: os.path.dirname(output.tab)
    conda:
        "../envs/krona.yml"
    singularity:
        "docker://biocontainers/krona:v2.7.1_cv1"
    shell:
        """
        ktUpdateTaxonomy.sh {params.taxdir} >{log} 2>&1
        """

rule classifier2krona:
    input:
        opj(config["paths"]["results"], "{classifier}",
            "{sample}_{unit}_{seq_type}.kreport"),
        opj("resources", "krona", "taxonomy.tab")
    output:
        opj(config["paths"]["results"], "{classifier}",
            "{sample}_{unit}_{seq_type}.html")
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
            {input[0]},{wildcards.sample}_{wildcards.unit}
        """

rule all2krona:
    input:
        f=get_all_files(samples, opj(config["paths"]["results"],
                                    "{classifier}"), ".kreport"),
        h=get_all_files(samples, opj(config["paths"]["results"],
                                    "{classifier}"), ".html"),
        t=opj("resources", "krona", "taxonomy.tab")
    output:
        opj(config["paths"]["results"], "report", "{classifier}", "{classifier}.krona.html")
    log:
        opj(config["paths"]["results"], "report", "{classifier}", "{classifier}.krona.log")
    params:
        tax="resources/krona",
        input_string=krona_input(config, samples, "{classifier}")
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
