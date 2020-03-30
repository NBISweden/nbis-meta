localrules:
    download_kraken_build

def get_kraken_index_url(config):
    url_base="ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/"
    d ={'minikraken2_v1_8GB': 'minikraken2_v1_8GB_201904_UPDATE.tgz',
         'minikraken2_v2_8GB': 'minikraken2_v2_8GB_201904_UPDATE.tgz',
         '16S_Greengenes': '16S_Greengenes_20190418.tgz',
         '16S_RDP': '16S_RDP_20190418.tgz',
         '16S_Silva': '16S_Silva_20190418.tgz'}
    try:
        url="{}/{}".format(url_base, d[config["kraken_prebuilt"]])
    except KeyError:
        url=""
    return url

rule download_kraken_build:
    """Downloads pre-built kraken2 index"""
    output:
        expand(opj(config["kraken_index_path"],"{n}.k2d"),
               n=["hash", "opts", "taxo"])
    params:
        dir=config["kraken_index_path"],
        tar=opj(config["scratch_path"],
                "{base}.tgz".format(base=config["kraken_prebuilt"])),
        url=get_kraken_index_url(config),
        tmpdir=opj(config["scratch_path"],"kraken_db")
    shell:
         """
         mkdir -p {params.tmpdir}
         curl -L -s -o {params.tar} {params.url}
         tar -C {params.tmpdir} -xvf {params.tar}
         mv {params.tmpdir}/*/* {params.dir}
         rm -r {params.tar} {params.tmpdir}  
         """