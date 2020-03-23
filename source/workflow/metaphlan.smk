# Workflow module for the Metaphlan2 classifier
metaphlan_input = get_all_files(samples,
                                opj(config["results_path"],"metaphlan2"),
                                ".mp2.out")
metaphlan_input.append(opj(config["report_path"],
                           "metaphlan2","metaphlan2.krona.html"))
metaphlan_input.append(opj(config["report_path"],
                           "metaphlan2","metaphlan2.heatmap.png"))

include: "../rules/Classify/metaphlan2.smk"