preprocess_input = [opj(config["report_path"],"samples_report.html")]

include: "../rules/Preprocess/preprocessing.smk"
include: "../rules/Preprocess/sample_report.smk"