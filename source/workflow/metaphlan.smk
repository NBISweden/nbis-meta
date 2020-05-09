# Workflow module for the Metaphlan2 classifier
metaphlan_input = get_all_files(samples,
                                opj(config["results_path"],"metaphlan"),
                                ".profile", nested=True)

include: "../rules/Classify/metaphlan.smk"