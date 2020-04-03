localrules:
    download_humann_uniref,
    download_chocophlan

rule build_metaphlan:
    """
    Download and build the metaphlan bowtie2 database
    """
    output:
        expand(opj(config["resource_path"], "metaphlan", "mpa_v20_m200.{s}.bt2"),
                   s = ["1","2","3","4","rev.1","rev.2"])
    log:
        opj(config["resource_path"], "metaphlan", "mpa.log")
    params:
        dir = opj(config["resource_path"], "metaphlan")
    threads: 4
    conda:
        "../../../envs/metaphlan.yml"
    shell:
        """
        metaphlan2.py --install --bowtie2db {params.dir} \
            --nproc {threads} >{log} 2>&1
        """

rule download_chocophlan:
    """
    Use the humann2_databases tool to download the required chocophlan database
    """
    output:
        touch(opj(config["resource_path"], "humann2", "chocophlan", "done"))
    log:
        opj(config["resource_path"], "humann2", "chocophlan", "download.log")
    params:
        dir = opj(config["resource_path"], "humann2", "chocophlan")
    conda:
        "../../../envs/metaphlan.yml"
    shell:
        """
        humann2_databases --download chocophlan full {params.dir} > {log} 2>&1
        """

rule download_humann_uniref:
    output:
        opj(config["resource_path"], "humann2", "uniref",
            "uniref{i}_annotated.1.1.dmnd")
    log:
        opj(config["resource_path"], "humann2", "uniref", "uniref{i}.log")
    params:
        dir = opj(config["resource_path"], "humann2")
    conda:
        "../../../envs/metaphlan.yml"
    shell:
        """
        humann2_databases \
            --download uniref uniref{wildcards.i}_diamond \
            {params.dir} >{log} 2>&1
        """