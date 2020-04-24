localrules:
    download_rfams,
    press_rfams

rule download_rfams:
    output:
        tar = temp(opj(config["infernal_dbpath"],"Rfam.tar.gz")),
        cm = opj(config["infernal_dbpath"],"Rfam.rRNA.cm"),
        readme = opj(config["infernal_dbpath"],"README"),
        version = opj(config["infernal_dbpath"],"Rfam.version"),
        clanin = opj(config["infernal_dbpath"],"Rfam.clanin")
    log:
        opj(config["infernal_dbpath"], "download.log")
    params:
        url = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT",
        rfams = " ".join(["RF00001.cm","RF00002.cm","RF00177.cm",
               "RF01959.cm","RF01960.cm","RF02540.cm","RF02541.cm",
               "RF02542.cm","RF02543.cm","RF02545.cm","RF02546.cm",
               "RF02547.cm","RF02554.cm","RF02555.cm"]),
        dir = lambda w, output: os.path.dirname(output.cm)
    shell:
        """
        # Get the Rfam tarball
        curl -o {output.tar} {params.url}/Rfam.tar.gz > {log} 2>&1
        # Extract only rfams of interest
        tar -C {params.dir} -zxf {output.tar} {params.rfams}
        cat {params.dir}/*.cm > {output.cm}
        
        # Get release
        curl -o {output.readme} {params.url}/README 2>/dev/null
        grep -m 1 Release {output.readme} > {output.version}
        
        # Get clans
        curl {params.url}/Rfam.clanin 2>/dev/null | egrep -w \
            "CL0011[123]" > {output.clanin}
        """

rule press_rfams:
    input:
        opj(config["infernal_dbpath"],"Rfam.rRNA.cm")
    output:
        expand(opj(config["infernal_dbpath"],"Rfam.rRNA.cm.i1{suffix}"),
               suffix=["m","i","f","p"])
    log:
        opj(config["infernal_dbpath"], "cmpress.log")
    conda:
        "../../../envs/annotation.yml"
    shell:
        """
        cmpress {input} > {log} 2>&1
        """