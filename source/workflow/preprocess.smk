preprocess_input = [opj(config["report_path"],"samples_report.html")]

include: "../rules/Preprocess/preprocessing.rules"
include: "../rules/Preprocess/sample_report.rules"