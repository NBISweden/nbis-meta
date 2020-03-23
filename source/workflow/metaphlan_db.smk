# Workflow module for generating the Metaphlan2 database

metaphlan_db_input=expand(opj(config["resource_path"],"metaphlan2",
                              "mpa_v20_m200.{index}.bt2"), index = range(1,5))

include: "../rules/Classify/prepare_metaphlan.smk"