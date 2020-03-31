# Workflow module for kraken2 classifier

kraken_input = get_all_files(samples, opj(config["results_path"],"kraken"),
                             ".kreport")
kraken_input.append(opj(config["report_path"],"kraken","kraken.krona.html"))

include: "../rules/Classify/kraken.smk"