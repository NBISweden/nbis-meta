localrules:
    metaphlan2kronatext,
    metaphlan2krona,
    all_metaphlan2_to_krona,
    merge_metaphlan2,

rule metaphlan_pe:
    input:
        R1 = opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R1"+PREPROCESS+".fastq.gz"),
        R2 = opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R2"+PREPROCESS+".fastq.gz"),
        db = expand(opj(config["resource_path"], "metaphlan",
                                "mpa_v20_m200.{s}.bt2"),
                   s = ["1","2","3","4","rev.1","rev.2"])
    output:
        opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_pe.profile")
    log:
        opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_pe.log")
    params:
        bt2 = opj(config["results_path"],"metaphlan","{sample}_{run}","bt2"),
        dir = opj(config["resource_path"], "metaphlan")
    conda:
        "../../../envs/metaphlan.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metaphlan2.py \
            {input.R1},{input.R2} \
            --bowtie2db {params.dir} \
            --bt2_ps very-sensitive \
            -t rel_ab \
            --bowtie2out {params.bt2} \
            --no_map \
            --nproc {threads} \
            --input_type multifastq \
            -o {output[0]} > {log} 2>&1
        """

rule metaphlan_se:
    input:
        se=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_se"+PREPROCESS+".fastq.gz")
    output:
        opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_se.profile")
    log:
        opj(config["results_path"],"metaphlan","{sample}_{run}",
            "{sample}_{run}_se.log")
    conda:
        "../../../envs/metaphlan.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metaphlan2.py \
            {input.se} \
            --bt2_ps very-sensitive \
            -t rel_ab \
            --nproc {threads} \
            --input_type fastq \
            -o {output[0]} >{log} 2>&1
        """

########################
## Create Krona plots ##
########################
rule metaphlan2kronatext:
    input:
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_{seq_type}.mp2.out.krona")
    output:
        temp(opj(config["results_path"],"metaphlan2",
                 "{sample}_{run}_{seq_type}.mp2.krona"))
    conda: 
        "../../../envs/metaphlan.yml"
    shell:
        """
        metaphlan2krona.py -p {input[0]} -k {output[0]}
        """

rule metaphlan2krona:
    input:
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_{seq_type}.mp2.krona")
    output:
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_{seq_type}.mp2.html")
    conda:
        "../../../envs/krona.yml"
    shell:
        """
        ktImportText \
            {input[0]},{wildcards.sample}_{wildcards.run} -o {output[0]}
        """

def get_mp_input_string(samples):
    input_string=""
    files=get_all_files(samples,opj(config["results_path"],
                                    "metaphlan2"),".mp2.krona")
    for f in files:
        sample_run=os.path.basename(f).replace("_pe.mp2.krona","").replace("_se.mp2.krona","")
        input_string+=" {},{}".format(f,sample_run)
    return input_string

rule all_metaphlan2_to_krona:
    """Combined krona plot for all samples"""
    input:
        f=get_all_files(samples, opj(config["results_path"],"metaphlan2"),".mp2.krona"),
        h=get_all_files(samples, opj(config["results_path"],"metaphlan2"),".mp2.html")
    output:
        opj(config["report_path"], "metaphlan2", "metaphlan2.krona.html")
    params:
        input_string=get_mp_input_string(samples)
    conda:
        "../../../envs/krona.yml"
    shell:
        """
        ktImportText {params.input_string} -o {output[0]}
        """

####################
## Merge and plot ##
####################
rule merge_metaphlan2:
    input:
        get_all_files(samples, opj(config["results_path"],"metaphlan"),
                                ".profile", nested=True)
    output:
        opj(config["report_path"],"metaphlan","metaphlan.merged.tsv")
    conda: 
        "../../../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} | sed 's/_[ps]e.mp2//g' > {output[0]}
        """