# Workflow module for kraken2 classifier

kraken_input = get_all_files(samples, opj(config["results_path"],"kraken"),
                             ".kreport")
# Exclude krona output until krona conda env is fixed on macos
if platform.uname().system != "Darwin":
    kraken_input.append(opj(config["report_path"],"kraken","kraken.krona.html"))

include: "../rules/Classify/kraken.smk"