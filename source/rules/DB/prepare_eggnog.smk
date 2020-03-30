localrules:
    download_eggnog,
    get_kegg_info

rule download_eggnog:
    output:
        opj(config["resource_path"],"eggnog-mapper","eggnog.db"),
        opj(config["resource_path"],"eggnog-mapper","eggnog.version")
    log:
        opj(config["resource_path"],"eggnog-mapper","download.log")
    conda:
        "../../../envs/annotation.yml"
    params:
        dbs="none",
        resource_path=opj(config["resource_path"], "eggnog-mapper")
    shell:
        """
        download_eggnog_data.py --data_dir \
            {params.resource_path} -y > {log} 2>&1
        egrep -o "emapperdb-[0-9].[0-9].[0-9]" {log} > {output[1]}
        """

rule get_kegg_info:
    #TODO: Check which files are needed with new eggnog-mapper version
    output:
        expand(opj(config["resource_path"],"kegg", "{f}"),
        f=["kegg_ec2pathways.tsv", 
             "kegg_ko2ec.tsv",  
             "kegg_ko2pathways.tsv", 
             "kegg_kos.tsv",
             "kegg_modules.tsv", 
             "kegg_pathways.tsv"])
    params:
        resource_dir=opj(config["resource_path"], "kegg")
    shell:
        """
        python source/utils/eggnog-parser.py download {params.resource_dir}
        """