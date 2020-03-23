# Workflow module for centrifuge classifier
centrifuge_input = get_all_files(samples,
                                 opj(config["results_path"],"centrifuge"),
                                ".kreport")
centrifuge_input.append(opj(config["report_path"],"centrifuge",
                            "centrifuge.krona.html"))

include: "../rules/Classify/centrifuge.smk"
