from scripts.common import classify_input, krona_input, metaphlan_krona_string

localrules:
    classify,
    krona_taxonomy,
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
        expand("resources/kraken/prebuilt/{kraken_prebuilt}/{n}.k2d",
               kraken_prebuilt=config["kraken"]["prebuilt"], n=["hash", "opts", "taxo"])
    log:
        "resources/kraken/prebuilt/{kraken_prebuilt}/download.log".format(kraken_prebuilt=config["kraken"]["prebuilt"])
    params:
        dir=lambda w, output: os.path.dirname(output[0]),
        tar="{temp}/{base}.tgz".format(temp=config["paths"]["temp"],
                                       base=config["kraken"]["prebuilt"]),
        url=config["kraken"]["prebuilt_url"],
        tmpdir = "{}/kraken_db".format(config["paths"]["temp"])
    shell:
         """
         mkdir -p {params.tmpdir}
         curl -L -o {params.tar} {params.url} > {log} 2>&1
         tar -C {params.tmpdir} -xf {params.tar}
         hashfile=$(find {params.tmpdir}/ -name "hash.k2d")
         hashdir=$(dirname $hashfile)
         mv $hashdir/* {params.dir}
         rm -r {params.tar} {params.tmpdir}
         """

rule kraken_build_standard:
    output:
        expand("resources/kraken/standard/{n}.k2d",
               n=["hash", "opts", "taxo"])
    log:
        build="resources/kraken/standard/build.log",
        clean="resources/kraken/standard/clean.log",
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
        R1=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_R1{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        R2=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_R2{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        db=expand("{kraken_index_path}/{n}.k2d",
                  kraken_index_path=config["kraken"]["index_path"],
                  n=["hash", "opts", "taxo"])
    output:
        expand("{results_path}/kraken/{{sample}}_{{unit}}_pe.out",
               results_path=config["paths"]["results"]),
        expand("{results_path}/kraken/{{sample}}_{{unit}}_pe.kreport",
               results_path=config["paths"]["results"])
    log:
        expand("{results_path}/kraken/{{sample}}_{{unit}}_pe.log",
               results_path=config["paths"]["results"])
    params:
        db=config["kraken"]["index_path"],
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
        se=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_se{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        db=expand("{kraken_index_path}/{n}.k2d",
                  kraken_index_path=config["kraken"]["index_path"],
                  n=["hash", "opts", "taxo"])
    output:
        expand("{results_path}/kraken/{{sample}}_{{unit}}_se.out",
               results_path=config["paths"]["results"]),
        expand("{results_path}/kraken/{{sample}}_{{unit}}_se.kreport",
               results_path=config["paths"]["results"])
    log:
        expand("{results_path}/kraken/{{sample}}_{{unit}}_se.log",
               results_path=config["paths"]["results"])
    params:
        db=config["kraken"]["index_path"],
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

rule kraken_contigs:
    input:
        expand("{results_path}/assembly/{{assembly}}/final_contigs.fa",
            results_path = config["paths"]["results"])
    output:
        expand("{results_path}/annotation/{{assembly}}/final_contigs.kraken.out",
            results_path=config["paths"]["results"]),
        expand("{results_path}/annotation/{{assembly}}/final_contigs.kraken.kreport",
           results_path=config["paths"]["results"])
    log:
        expand("{results_path}/annotation/{{assembly}}.kraken.log",
            results_path=config["paths"]["results"])
    params:
        db=config["kraken"]["index_path"],
        mem=config["kraken"]["mem"]
    threads: 10
    resources:
        runtime= lambda wildcards,attempt: attempt ** 2 * 60 * 10
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
        db=expand("{centrifuge_dir}/{base}.{i}.cf",
                  centrifuge_dir=config["centrifuge"]["dir"],
                  base=config["centrifuge"]["base"],
                  i=[1, 2, 3],)
    log:
        "{centrifuge_dir}/download.log".format(centrifuge_dir=config["centrifuge"]["dir"])
    params:
        dir=config["centrifuge"]["dir"],
        tar="{centrifuge_dir}/{base}.tar.gz".format(base=config["centrifuge"]["base"],
                                                    centrifuge_dir=config["centrifuge"]["dir"]),
        url=config["centrifuge"]["prebuilt_url"]
    shell:
        """
        curl -o {params.tar} {params.url} > {log} 2>&1
        tar -C {params.dir} -xf {params.tar} >>{log} 2>&1
        rm {params.tar}
        """

rule centrifuge_pe:
    input:
        R1=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_R1{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        R2=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_R2{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        db=expand("{centrifuge_dir}/{base}.{i}.cf",
                  centrifuge_dir=config["centrifuge"]["dir"],
                  base=config["centrifuge"]["base"], i=[1, 2, 3])
    output:
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_pe.out",
               results_path=config["paths"]["results"]),
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_pe.report",
               results_path=config["paths"]["results"])
    log:
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_pe.log",
               results_path=config["paths"]["results"])
    params:
        prefix="{centrifuge_dir}/{base}".format(centrifuge_dir=config["centrifuge"]["dir"],
                                                base=config["centrifuge"]["base"]),
        tmp_out=expand("{temp}/{{sample}}_{{unit}}_pe.out",
                       temp=config["paths"]["temp"]),
        tmpdir=config["paths"]["temp"],
        tmp_report=expand("{temp}/{{sample}}_{{unit}}_pe.report",
                          temp=config["paths"]["temp"]),
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
        se=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_se{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        db=expand("{centrifuge_dir}/{base}.{i}.cf",
                  centrifuge_dir=config["centrifuge"]["dir"],
                  base=config["centrifuge"]["base"], i=[1, 2, 3])
    output:
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_se.out",
               results_path=config["paths"]["results"]),
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_se.report",
               results_path=config["paths"]["results"])
    log:
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_se.log",
               results_path=config["paths"]["results"])
    params:
        prefix="{centrifuge_dir}/{base}".format(centrifuge_dir=config["centrifuge"]["dir"],
                                                base=config["centrifuge"]["base"]),
        tmp_out=expand("{temp}/{{sample}}_{{unit}}_se.out",
                       temp=config["paths"]["temp"]),
        tmpdir=config["paths"]["temp"],
        tmp_report=expand("{temp}/{{sample}}_{{unit}}_se.report",
                          temp=config["paths"]["temp"]),
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
        f=expand("{results_path}/centrifuge/{{sample}}_{{unit}}_{{seq_type}}.out",
            results_path=config["paths"]["results"]),
        db=expand("{dirname}/{base}.{i}.cf",
                  dirname=config["centrifuge"]["dir"],
                  i=[1, 2, 3], base=config["centrifuge"]["base"])
    output:
        expand("{results_path}/centrifuge/{{sample}}_{{unit}}_{{seq_type}}.kreport",
            results_path=config["paths"]["results"])
    params:
        min_score=config["centrifuge"]["min_score"],
        prefix="{dirname}/{base}".format(dirname=config["centrifuge"]["dir"],
                                         base=config["centrifuge"]["base"])
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
        expand("resources/metaphlan/{index}.{s}.bt2",
               index=config["metaphlan"]["index"],
               s=["1", "2", "3", "4", "rev.1", "rev.2"])
    log:
        "resources/metaphlan/mpa.log"
    params:
        dir="resources/metaphlan",
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
        R1=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_R1{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        R2=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_R2{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        db=expand("resources/metaphlan/{index}.{s}.bt2",
                  index=config["metaphlan"]["index"],
                  s=["1", "2", "3", "4", "rev.1", "rev.2"])
    output:
        tsv=expand("{results_path}/metaphlan/{{sample}}_{{unit}}_pe.tsv",
                   results_path=config["paths"]["results"]),
        bt2=expand("{results_path}/metaphlan/{{sample}}_{{unit}}_pe.bt2",
                   results_path=config["paths"]["results"])
    log:
        expand("{results_path}/metaphlan/{{sample}}_{{unit}}_pe.log",
               results_path=config["paths"]["results"])
    params:
        dir="resources/metaphlan"
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
        se=expand("{results_path}/intermediate/preprocess/{{sample}}_{{unit}}_se{preprocess}.fastq.gz",
                  results_path=config["paths"]["results"], preprocess=PREPROCESS),
        db=expand("resources/metaphlan/{index}.{s}.bt2",
                  index=config["metaphlan"]["index"],
                  s=["1", "2", "3", "4", "rev.1", "rev.2"])
    output:
        tsv=expand("{results_path}/metaphlan/{{sample}}_{{unit}}_se.tsv",
                   results_path=config["paths"]["results"]),
        bt2=expand("{results_path}/metaphlan/{{sample}}_{{unit}}_se.bt2",
                   results_path=config["paths"]["results"])
    log:
        expand("{results_path}/metaphlan/{{sample}}_{{unit}}_se.log",
               results_path=config["paths"]["results"])
    params:
        dir = "resources/metaphlan"
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
        get_all_files(samples, "{}/metaphlan".format(config["paths"]["results"]), ".tsv")
    output:
        "{}/report/metaphlan/metaphlan.tsv".format(config["paths"]["results"])
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """

rule metaphlan2krona_table:
    input:
        expand("{results_path}/metaphlan/{{sample}}_{{unit}}_{{seq_type}}.tsv",
               results_path=config["paths"]["results"])
    output:
        temp(expand("{results_path}/metaphlan/{{sample}}_{{unit}}_{{seq_type}}.krona",
                    results_path=config["paths"]["results"]))
    script:
        "../scripts/classification_utils.py"

rule metaphlan2krona:
    input:
        files=get_all_files(samples, "{}/metaphlan".format(config["paths"]["results"]), ".krona"),
        db="resources/krona/taxonomy.tab"
    output:
        "{}/report/metaphlan/metaphlan.html".format(config["paths"]["results"])
    log:
        "{}/report/metaphlan/krona.log".format(config["paths"]["results"])
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
        "{}/report/metaphlan/metaphlan.tsv".format(config["paths"]["results"])
    output:
        "{}/report/metaphlan/metaphlan.pdf".format(config["paths"]["results"])
    params:
        rank=config["metaphlan"]["plot_rank"]
    conda:
        "../envs/plotting.yml"
    notebook:
        "../notebooks/metaphlan.py.ipynb"

##### krona #####

rule krona_taxonomy:
    output:
        tab="resources/krona/taxonomy.tab"
    log:
        "resources/krona/taxonomy.log"
    params:
        taxdir=lambda w, output: os.path.dirname(output.tab)
    conda:
        "../envs/krona.yml"
    container:
        "docker://biocontainers/krona:v2.7.1_cv1"
    shell:
        """
        ktUpdateTaxonomy.sh {params.taxdir} >{log} 2>&1
        """

rule classifier2krona:
    input:
        expand("{results_path}/{{classifier}}/{{sample}}_{{unit}}_{{seq_type}}.kreport",
               results_path = config["paths"]["results"]),
        "resources/krona/taxonomy.tab"
    output:
        expand("{results_path}/{{classifier}}/{{sample}}_{{unit}}_{{seq_type}}.html",
               results_path = config["paths"]["results"])
    log:
        expand("{results_path}/{{classifier}}/{{sample}}_{{unit}}_{{seq_type}}.krona.log",
               results_path = config["paths"]["results"])
    params:
        tax="resources/krona"
    conda:
        "../envs/krona.yml"
    container:
        "docker://continuumio/miniconda3:4.8.2"
    shell:
        """
        ktImportTaxonomy -t 5 -m 3 -tax {params.tax} -o {output[0]} \
            {input[0]},{wildcards.sample}_{wildcards.unit} > {log} 2>&1
        """

rule all2krona:
    input:
        f=get_all_files(samples,
                        expand("{results_path}/{{classifier}}", results_path = config["paths"]["results"]),
                        ".kreport"),
        h=get_all_files(samples,
                        expand("{results_path}/{{classifier}}", results_path = config["paths"]["results"]),
                        ".html"),
        t="resources/krona/taxonomy.tab"
    output:
        expand("{results_path}/report/{{classifier}}/{{classifier}}.krona.html",
               results_path = config["paths"]["results"])
    log:
        expand("{results_path}/report/{{classifier}}/{{classifier}}.krona.log",
               results_path = config["paths"]["results"])
    params:
        tax="resources/krona",
        input_string=krona_input(config, samples, "{classifier}")
    conda:
        "../envs/krona.yml"
    container:
        "docker://continuumio/miniconda3:4.8.2"
    shell:
         """
         ktImportTaxonomy \
            -t 5 -m 3 -tax {params.tax} -o {output[0]} \
            {params.input_string} > {log} 2>&1
         """
