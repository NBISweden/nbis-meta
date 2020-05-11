# Workflow module for generating the Metaphlan2 database

metaphlan_db_input = expand(opj(config["resource_path"], "metaphlan", "{index}.{s}.bt2"),
                            index = config["metaphlan_index"],
                            s = ["1","2","3","4","rev.1","rev.2"])

include: "../rules/Classify/prepare_metaphlan.smk"