# Workflow module for generating the Metaphlan2 database

metaphlan_db_input = opj(config["resource_path"],"metaphlan","bowtie2-build.log")

include: "../rules/Classify/prepare_metaphlan.smk"