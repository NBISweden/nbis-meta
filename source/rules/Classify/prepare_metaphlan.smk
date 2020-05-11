rule build_metaphlan:
    """
    Download and build the metaphlan bowtie2 database
    """
    output:
        expand(opj(config["resource_path"], "metaphlan", "{index}.{s}.bt2"),
               index = config["metaphlan_index"], s = ["1","2","3","4","rev.1","rev.2"])
    log:
        opj(config["resource_path"], "metaphlan", "mpa.log")
    params:
        dir = opj(config["resource_path"], "metaphlan"),
        index = config["metaphlan_index"]
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*1
    conda:
        "../../../envs/metaphlan.yml"
    shell:
        """
        metaphlan --install --bowtie2db {params.dir} \
            --nproc {threads} -x {params.index} >{log} 2>&1
        """