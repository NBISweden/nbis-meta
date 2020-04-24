localrules:
    download_eggnog,
    get_kegg_info

rule download_eggnog:
    output:
        db = opj(config["resource_path"],"eggnog-mapper","eggnog.db"),
        version = opj(config["resource_path"],"eggnog-mapper","eggnog.version")
    log:
        opj(config["resource_path"],"eggnog-mapper","download.log")
    conda:
        "../../../envs/annotation.yml"
    params:
        dbs="none",
        data_dir = lambda w, output: os.path.dirname(output.db)
    shell:
        """
        download_eggnog_data.py --data_dir {params.data_dir} -y > {log} 2>&1
        egrep -o "emapperdb-[0-9].[0-9].[0-9]" {log} > {output.version}
        """

rule get_kegg_info:
    #TODO: Check which files are needed with new eggnog-mapper version
    output:
        expand(opj(config["resource_path"], "kegg", "{f}"),
        f=["kegg_ec2pathways.tsv", 
             "kegg_ko2ec.tsv",  
             "kegg_ko2pathways.tsv", 
             "kegg_kos.tsv",
             "kegg_modules.tsv", 
             "kegg_pathways.tsv"])
    log:
        opj(config["resource_path"], "kegg", "download.log")
    params:
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        """
        python source/utils/eggnog-parser.py \
            download {params.outdir} > {log} 2>&1
        """