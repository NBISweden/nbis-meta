localrules:
    download_metaphlan

rule download_metaphlan:
    output:
        opj(config["resource_path"],"metaphlan2","mpa_v20_m200.fna"),
        opj(config["resource_path"],"metaphlan2","md5")
    params:
        outdir=opj(config["resource_path"],"metaphlan2"),
        tmpdir=opj(config["scratch_path"],"metaphlan2"),
        url_base="https://bitbucket.org/biobakery/metaphlan2/downloads"
    shell:
        """
        mkdir -p {params.tmpdir}
        curl -s -L -o {params.tmpdir}/mpa.tar {params.url_base}/mpa_v20_m200.tar
        tar -C {params.tmpdir} -xf {params.tmpdir}/mpa.tar
        rm {params.tmpdir}/mpa.tar 
        curl -s -L -o {params.tmpdir}/md5 {params.url_base}/mpa_v20_m200.md5
        bunzip2 {params.tmpdir}/mpa_v20_m200.fna.bz2
        mv {params.tmpdir}/* {params.outdir}
        """

rule bowtiebuild_metaphlan:
    input:
        opj(config["resource_path"],"metaphlan2","mpa_v20_m200.fna"),
        opj(config["resource_path"],"metaphlan2","md5")
    output:
        expand(opj(config["resource_path"],"metaphlan2",
                   "mpa_v20_m200.{index}.bt2"), index=range(1,5)),
        expand(opj(config["resource_path"],"metaphlan2",
                   "mpa_v20_m200.rev.{rev_index}.bt2"), rev_index=range(1,3))
    log:
        opj(config["resource_path"],"metaphlan2","bowtie_build.log")
    params:
        prefix=opj(config["resource_path"],"metaphlan2","mpa_v20_m200")
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*48
    conda:
        "../../../envs/metaphlan2.yml"
    shell:
        """
        bowtie2-build \
            --threads {threads} {input[0]} {params.prefix} >{log} 2>&1
        """

