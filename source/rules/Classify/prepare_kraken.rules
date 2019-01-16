localrules:
    download_kraken_build

def get_kraken_index_url(config):
    d = {'minikraken2_v1_8GB': 'http://ccb.jhu.edu/software/kraken2/dl/minikraken2_v1_8GB.tgz',
         'minikraken2_v2_8GB': 'http://ccb.jhu.edu/software/kraken2/dl/minikraken2_v2_8GB.tgz'}
    try:
        url = d[config["kraken_prebuilt"]]
    except KeyError:
        url = ""
    return url

rule download_kraken_build:
    """Downloads pre-built kraken2 index"""
    output:
        expand(opj(config["kraken_index_path"], "{n}.k2d"), n = ["hash", "opts", "taxo"])
    params:
        dir = config["kraken_index_path"],
        tar = opj(config["kraken_index_path"], "{base}.tgz".format(base=config["kraken_prebuilt"]))
    run:
        url = get_kraken_index_url(config)
        shell("curl -o {params.tar} {url}")
        shell("tar -C {params.dir} -xvf {params.tar}")
        shell("rm {params.tar}")