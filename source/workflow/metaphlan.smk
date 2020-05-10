# Workflow module for the Metaphlan2 classifier
metaphlan_input = [opj(config["report_path"], "metaphlan", "metaphlan.tsv")]

if platform.uname().system != "Darwin":
    metaphlan_input.append(opj(config["report_path"], "metaphlan",
                               "metaphlan.html"))

include: "../rules/Classify/metaphlan.smk"