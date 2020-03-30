localrules:
    centrifuge_kreport

rule centrifuge_pe:
    input:
        R1=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["centrifuge_dir"],"{base}.{i}.cf"),
                  i=[1,2,3], base=config["centrifuge_base"])
    output:
        opj(config["results_path"],"centrifuge","{sample}_{run}_pe.out"),
        opj(config["results_path"],"centrifuge","{sample}_{run}_pe.report")
    params:
        prefix=opj(config["centrifuge_dir"],
                   "{base}".format(base=config["centrifuge_base"])),
        tmp_out=opj(config["scratch_path"],"{sample}_{run}_pe.out"),
        tmp_report=opj(config["scratch_path"],"{sample}_{run}_pe.report")
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/centrifuge.yml"
    shell:
        """
        mkdir -p {config[scratch_path]}
        centrifuge \
            -k {config[centrifuge_max_assignments]} \
            -1 {input.R1} \
            -2 {input.R2} \
            -x {params.prefix} \
            -S {params.tmp_out} \
            --report-file {params.tmp_report} \
            -p {threads}
        mv {params.tmp_out} {output[0]}
        mv {params.tmp_report} {output[1]}
        """

rule centrifuge_se:
    input:
        se=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["centrifuge_dir"],
                      "{base}.{i}.cf"), i=[1,2,3],
                  base=config["centrifuge_base"])
    output:
        opj(config["results_path"],"centrifuge","{sample}_{run}_se.out"),
        opj(config["results_path"],"centrifuge","{sample}_{run}_se.report")
    params:
        prefix=opj(config["centrifuge_dir"],
                   "{base}".format(base=config["centrifuge_base"])),
        tmp_out=opj(config["scratch_path"],"{sample}_{run}_se.out"),
        tmp_report=opj(config["scratch_path"],"{sample}_{run}_se.report")
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/centrifuge.yml"
    shell:
        """
        mkdir -p {config[scratch_path]}
        centrifuge \
            -k {config[centrifuge_max_assignments]} \
            -U {input.se} \
            -x {params.prefix} \
            -S {params.tmp_out} \
            --report-file {params.tmp_report} \
            -p {threads}
        mv {params.tmp_out} {output[0]}
        mv {params.tmp_report} {output[1]}
        """

######################
## Generate reports ##
######################
rule centrifuge_kreport:
    input:
        f=opj(config["results_path"],"centrifuge",
              "{sample}_{run}_{seq_type}.out"),
        db=expand(opj(config["centrifuge_dir"],
                      "{base}.{i}.cf"), i=[1,2,3],
                  base=config["centrifuge_base"])
    output:
        opj(config["results_path"],"centrifuge",
            "{sample}_{run}_{seq_type}.kreport")
    params:
        min_score=config["centrifuge_min_score"],
        prefix=opj(config["centrifuge_dir"],
                   "{base}".format(base=config["centrifuge_base"]))
    conda:
        "../../../envs/centrifuge.yml"
    shell:
        """
        centrifuge-kreport \
            --min-score {params.min_score} \
            -x {params.prefix} \
            {input.f} > {output[0]}
        """