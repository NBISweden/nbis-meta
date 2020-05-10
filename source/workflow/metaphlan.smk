# Workflow module for the Metaphlan2 classifier
metaphlan_input = opj(config["report_path"], "metaphlan", "metaphlan.tsv")

include: "../rules/Classify/metaphlan.smk"