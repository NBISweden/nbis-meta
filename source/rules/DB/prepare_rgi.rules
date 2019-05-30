localrules: download_rgi_data

rule download_rgi_data:
    output:
        opj(config["resource_path"], "card", "card.json"),
        opj(config["resource_path"], "card", "card.version"),
    params:
        tar = opj(config["resource_path"], "card", "data.tar.gz"),
        dir = opj(config["resource_path"], "card")
    shell:
         """
         curl -L -o {params.tar} https://card.mcmaster.ca/latest/data
         tar -C {params.dir} -xvf {params.tar} ./card.json
         # Store download date in versionfile
         date > {output[1]}
         rm {params.tar}
         """
