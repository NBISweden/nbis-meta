localrules:
    samples_qc_report

rule fastqc:
    """Run fastqc on preprocessed data"""
    input:
        fastq = opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_{pair}"+PREPROCESS+".fastq.gz")
    output:
        zip = opj(config["intermediate_path"],"fastqc",
            "{sample}_{run}_{pair}"+PREPROCESS+"_fastqc.zip")
    log:
        opj(config["intermediate_path"],"fastqc",
            "{sample}_{run}_{pair}.log")
    params:
        dir = lambda w, output: os.path.dirname(output.zip)
    shadow: "shallow"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        fastqc -q --noextract -o {params.results_path} {input} >{log} 2>&1
        """

rule aggregate_logs:
    """Rule for aggregating preprocessing logs"""
    input:
        trimlogs=get_trim_logs,
        sortmernalogs=get_sortmerna_logs,
        filtlogs=get_filt_logs,
        fastqc=get_fastqc_files
    output:
        flag = touch(temp(opj(config["report_path"],"multiqc_input","flag")))
    log:
        opj(config["report_path"], ".aggregate.log")
    params:
        dir = lambda w, output: os.path.dirname(output.flag)
    shell:
        """
        cp {input} {params.dir} > {log} 2>&1
        """

rule samples_qc_report:
    """Summarize sample QC statistics in a report """
    input:
        flag = opj(config["report_path"],"multiqc_input","flag")
    output:
        html = opj(config["report_path"],"samples_report.html"),
        txt = opj(config["report_path"],"samples_report_data",
            "multiqc_general_stats.txt")
    log:
        opj(config["report_path"], "multiqc.log")
    shadow:
        "shallow"
    params:
        config="config/multiqc_preprocess_config.yaml",
        output_dir = lambda w, output: os.path.dirname(output.html),
        input_dir = lambda w, input: os.path.dirname(input.flag)
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        multiqc \
            -f \
            -c {params.config} \
            -n samples_report.html \
            -o {params.output_dir} \
            {params.input_dir} >{log} 2>{log}
        rm -r {params.input_dir}
        """