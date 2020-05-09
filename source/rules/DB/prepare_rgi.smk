localrules: download_rgi_data

rule download_rgi_data:
    output:
        json = opj(config["resource_path"], "card", "card.json"),
        version = opj(config["resource_path"], "card", "card.version")
    log:
        opj(config["resource_path"], "card", "log")
    params:
        tar = opj(config["resource_path"], "card", "data.tar.gz"),
        dir = lambda w, output: os.path.dirname(output.json)
    shell:
         """
         curl -L -o {params.tar} \
            https://card.mcmaster.ca/latest/data >{log} 2>&1
         tar -C {params.dir} -xf {params.tar} ./card.json
         # Store download date in versionfile
         date > {output.version}
         rm {params.tar}
         """
