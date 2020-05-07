localrules:
    download_kraken_build

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
        dir = lambda w, output: os.path.dirname(output[0]),
        tar = opj(config["scratch_path"],
                "{base}.tgz".format(base=config["kraken_prebuilt"])),
        url = get_kraken_index_url(config),
        tmpdir = opj(config["scratch_path"],"kraken_db")
    shell:
         """
         mkdir -p {params.tmpdir}
         curl -L -o {params.tar} {params.url} > {log} 2>&1
         tar -C {params.tmpdir} -xf {params.tar}
         mv {params.tmpdir}/*/* {params.dir}
         rm -r {params.tar} {params.tmpdir}  
         """

rule kraken_build_standard:
    output:
        expand(opj(config["resource_path"], "kraken", "standard", "{n}.k2d"),
               n=["hash", "opts", "taxo"])
    log:
        build = opj(config["resource_path"], "kraken", "standard", "build.log"),
        clean = opj(config["resource_path"], "kraken", "standard", "clean.log"),
    params:
        dir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../../../envs/kraken.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*24
    shell:
        """
        kraken2-build --standard --db {params.dir} --threads {threads} > {log.build} 2>&1
        kraken2-build --clean {params.dir} > {log.clean} 2>&1
        """