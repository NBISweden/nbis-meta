# Workflow module for generating the Metaphlan2 database

metaphlan_db_input = expand(opj(config["resource_path"], "metaphlan",
                                "mpa_v20_m200.{s}.bt2"),
                   s = ["1","2","3","4","rev.1","rev.2"])

include: "../rules/Classify/prepare_metaphlan.smk"