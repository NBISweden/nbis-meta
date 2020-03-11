localrules:
    download_rfams,
    press_rfams

rule download_rfams:
    output:
        temp(opj(config["infernal_dbpath"],"Rfam.tar.gz")),
        opj(config["infernal_dbpath"],"Rfam.rRNA.cm"),
        opj(config["infernal_dbpath"],"README"),
        opj(config["infernal_dbpath"],"Rfam.version"),
        opj(config["infernal_dbpath"],"Rfam.clanin")
    params:
        url="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT",
        rfams=" ".join(["RF00001.cm","RF00002.cm","RF00177.cm","RF01118.cm",
               "RF01959.cm","RF01960.cm","RF02540.cm","RF02541.cm",
               "RF02542.cm","RF02543.cm","RF02545.cm","RF02546.cm",
               "RF02547.cm","RF02554.cm","RF02555.cm"]),
        dir=opj(config["infernal_dbpath"])
    shell:
        """
        # Get the Rfam tarball
        curl -o {output[0]} {params.url}/Rfam.tar.gz
        # Extract only rfams of interest
        tar -C {params.dir} -zxf {output[0]} {params.rfams}
        cat {params.dir}/*.cm > {output[1]}
        
        # Get release
        curl -o {output[2]} {params.url}/README
        grep -m 1 Release {output[2]} > {output[3]}
        
        # Get clans
        curl {params.url}/Rfam.clanin | egrep -w "CL0011[123]" > {output[4]}
        """

rule press_rfams:
    input:
        opj(config["infernal_dbpath"],"Rfam.rRNA.cm")
    output:
        expand(opj(config["infernal_dbpath"],"Rfam.rRNA.cm.i1{suffix}"),
               suffix=["m","i","f","p"])
    conda:
        "../../../envs/annotation.yml"
    shell:
        """
        cmpress {input}
        """