localrules:
    samples_qc_report

rule fastqc:
    """Run fastqc on preprocessed data"""
    input:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_{pair}"+PREPROCESS+".fastq.gz")
    output:
        opj(config["intermediate_path"],"fastqc",
            "{sample}_{run}_{pair}"+PREPROCESS+"_fastqc.zip")
    params:
        results_path=opj(config["intermediate_path"],"fastqc")
    shadow: "shallow"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/preprocess.yml"
    shell:
        "fastqc -q --noextract -o {params.results_path} {input}"

rule aggregate_logs:
    """Rule for aggregating preprocessing logs"""
    input:
        trimlogs=get_trim_logs,
        sortmernalogs=get_sortmerna_logs,
        filtlogs=get_filt_logs,
        fastqc=get_fastqc_files
    output:
        touch(temp(opj(config["report_path"],"multiqc_input","flag")))
    params:
        output_dir=opj(config["report_path"],"multiqc_input")
    run:
        for file in input.trimlogs:
            shell("cp {file} {params.output_dir}")
        for file in input.sortmernalogs:
            shell("cp {file} {params.output_dir}")
        for file in input.filtlogs:
            shell("cp {file} {params.output_dir}")
        for file in input.fastqc:
            shell("cp {file} {params.output_dir}")

rule samples_qc_report:
    """Summarize sample QC statistics in a report """
    input:
        opj(config["report_path"],"multiqc_input","flag")
    output:
        opj(config["report_path"],"samples_report.html"),
        opj(config["report_path"],"samples_report_data",
            "multiqc_general_stats.txt")
    log: opj(config["report_path"], "multiqc.log")
    shadow:
        "shallow"
    params:
        config="config/multiqc_preprocess_config.yaml",
        output_dir=opj(config["report_path"]),
        input_dir=opj(config["report_path"],"multiqc_input")
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        multiqc \
            -f \
            -c {params.config} \
            -n samples_report.html \
            -o {params.output_dir} \
            $(dirname {input}) >{log} 2>{log}
        rm -r {params.input_dir}
        """