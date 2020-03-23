# Workflow module for annotation database files

db_input = []

if config["taxonomic_annotation"]:
    db_input.append(opj(config["resource_path"],config["taxdb"],"diamond.dmnd".format(config["taxdb"])))
    db_input.append(opj(config["resource_path"], "taxonomy/taxonomy.sqlite"))
if config["infernal"]:
    db_input+=expand(opj(config["infernal_dbpath"],"Rfam.rRNA.cm.i1{suffix}"),
               suffix=["m","i","f","p"])
if config["eggnog"]:
    db_input.append(opj(config["resource_path"],"eggnog-mapper","eggnog.db"))
    db_input.append(opj(config["resource_path"],"eggnog-mapper","eggnog.version"))
    db_input += expand(opj(config["resource_path"],"kegg", "{f}"),
        f = ["kegg_ec2pathways.tsv", "kegg_ko2ec.tsv", "kegg_ko2modules.tsv", "kegg_ko2pathways.tsv", "kegg_kos.tsv",
             "kegg_modules.tsv", "kegg_pathways.tsv"])
if config["pfam"]:
    db_input.append(opj(config["resource_path"],"pfam","Pfam-A.hmm"))
    db_input.append(opj(config["resource_path"],"pfam","Pfam-A.hmm.h3f"))
    db_input.append(opj(config["resource_path"],"pfam","Pfam-A.hmm.dat"))
    db_input.append(opj(config["resource_path"],"pfam","Pfam-A.clans.tsv"))
if config["sortmerna"]:
    for f in config["sortmerna_dbs"]:
        db_input.append(opj(config["resource_path"],"rRNA_databases","{}.idx.stats".format(f)))
if config["rgi"]:
    db_input.append(opj(config["resource_path"], "card", "card.json"))

include: "../rules/DB/prepare_sortmerna.smk"
include: "../rules/DB/prepare_eggnog.rules"
include: "../rules/DB/prepare_hmms.smk"
include: "../rules/DB/prepare_infernal.smk"
include: "../rules/DB/prepare_taxonomy.smk"
include: "../rules/DB/prepare_rgi.rules"