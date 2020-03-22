localrules:
    download_centrifuge_build

def get_centrifuge_index_url(config):
    url_base="ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data"
    d={'nt_2018_2_12': 'nt_2018_2_12.tar.gz',
         'nt_2018_3_3': 'nt_2018_3_3.tar.gz',
         'p+h+v': 'p+h+v.tar.gz',
         'p_compressed+h+v': 'p_compressed+h+v.tar.gz',
         'p_compressed_2018_4_15': 'p_compressed_2018_4_15.tar.gz'}
    try:
        url="{}/{}".format(url_base, d[config["centrifuge_base"]])
    except KeyError:
        url=""
    return url

rule download_centrifuge_build:
    """Downloads pre-built centrifuge index"""
    output:
        db=expand(opj(config["centrifuge_dir"],
                      "{base}.{i}.cf"), i=[1,2,3],
                  base=config['centrifuge_base'])
    params:
        dir=config["centrifuge_dir"],
        tar=opj(config["centrifuge_dir"],
                "{base}.tar.gz".format(base=config["centrifuge_base"])),
        url=get_centrifuge_index_url(config)
    shell:
        """
        curl -s -o {params.tar} {url}
        tar -C {params.dir} -xvf {params.tar}
        rm {params.tar}
        """