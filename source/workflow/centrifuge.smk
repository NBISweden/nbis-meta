centrifuge_input = []

centrifuge_input.append(opj(config["report_path"],"centrifuge",
                            "centrifuge.krona.html"))

include: "../rules/Classify/centrifuge.smk"
