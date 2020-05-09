localrules:
    download_centrifuge_build

rule download_centrifuge_build:
    """Downloads pre-built centrifuge index"""
    output:
        db=expand(opj(config["centrifuge_dir"],
                      "{base}.{i}.cf"), i=[1,2,3],
                  base=config['centrifuge_base'])
    log:
        opj(config["centrifuge_dir"], "download.log")
    params:
        dir=config["centrifuge_dir"],
        tar=opj(config["centrifuge_dir"],
                "{base}.tar.gz".format(base=config["centrifuge_base"])),
        url=get_centrifuge_index_url(config)
    shell:
        """
        curl -o {params.tar} {params.url} > {log} 2>&1
        tar -C {params.dir} -xf {params.tar} >>{log} 2>&1
        rm {params.tar}
        """