localrules:
    metaphlan2kronatext,
    metaphlan2krona,
    all_metaphlan2_to_krona,
    merge_metaphlan2,
    metaphlan2_heatmap,
    metaphlan2graphlan,
    graphlan_annotate

rule metaphlan2_pe:
    input:
        R1=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R1"+PREPROCESS+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R2"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["resource_path"],"metaphlan2",
                        "mpa_v20_m200.{index}.bt2"), index=range(1,5))
    output:
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_pe.mp2.stats"),
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_pe.mp2.out"),
        temp(opj(config["results_path"],"metaphlan2",
                 "{sample}_{run}_pe.mp2.out.krona")),
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_pe.bowtie2.bz2")
    params:
        bowtie2db=opj(config["resource_path"],"metaphlan2")
    conda:
        "../../../envs/metaphlan2.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metaphlan2.py {input.R1},{input.R2} --min_alignment_len 50 \
            --bt2_ps very-sensitive -t rel_ab_w_read_stats \
            --bowtie2db {params.bowtie2db} --bowtie2out {output[3]} \
            --nproc {threads} \
            --input_type fastq > {output[0]}
        # Create the simple mp2 output
        grep "#SampleID" {output[0]} > {output[1]}
        egrep -v "clade_name|estimated" {output[0]} | cut -f1,2 >> {output[1]}
        # Create krona-compatible output
        grep "#SampleID" {output[0]} > {output[2]}
        grep "t__" {output[0]} | cut -f1,2 >> {output[2]}
        """

rule metaphlan2_se:
    input:
        se=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_se"+PREPROCESS+".fastq.gz"),
        db=expand(opj(config["resource_path"],"metaphlan2",
                        "mpa_v20_m200.{index}.bt2"), index=range(1,5))
    output:
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_se.mp2.stats"),
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_se.mp2.out"),
        temp(opj(config["results_path"],"metaphlan2",
                 "{sample}_{run}_se.mp2.out.krona")),
        opj(config["results_path"],"metaphlan2",
            "{sample}_{run}_se.bowtie2.bz2")
    params:
        bowtie2db=opj(config["resource_path"],"metaphlan2")
    conda: 
        "../../../envs/metaphlan2.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metaphlan2.py \
            {input.se} --min_alignment_len 50 \
            --bt2_ps very-sensitive --bowtie2db {params.bowtie2db} \
            -t rel_ab_w_read_stats --bowtie2out {output[3]} \
            --nproc {threads} --input_type fastq > {output[0]}
        # Create the simple mp2 output
        grep "#SampleID" {output[0]} > {output[1]}
        egrep -v "clade_name|estimated" {output[0]} | cut -f1,2 >> {output[1]}
        # Create krona-compatible output
        grep "#SampleID" {output[0]} > {output[2]}
        grep "t__" {output[0]} | cut -f1,2 >> {output[2]}
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
        "../../../envs/metaphlan2.yml"
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
        get_all_files(samples, opj(config["results_path"],"metaphlan2"),".mp2.out")
    output:
        opj(config["report_path"],"metaphlan2","metaphlan2.merged.out")
    conda: 
        "../../../envs/metaphlan2.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} | sed 's/_[ps]e.mp2//g' > {output[0]}
        """

rule metaphlan2_heatmap:
    #TODO: Due to deprecated matplotlib properties and incompatible conda
    #conda versions this rule needs to be updated with proper plotting
    input:
        opj(config["report_path"],"metaphlan2","metaphlan2.merged.out")
    output:
        touch(opj(config["report_path"],"metaphlan2","metaphlan2.heatmap.png"))
    shell:
        """
        """

rule metaphlan2graphlan:
    input:
        opj(config["report_path"],"metaphlan2","metaphlan2.merged.out")
    output:
        temp(opj(config["report_path"],"metaphlan2","metaphlan2.tree.txt")),
        temp(opj(config["report_path"],"metaphlan2","metaphlan2.annot.txt"))
    conda:
        "../../../envs/graphlan.yml"
    shell:
        """
        export2graphlan.py \
            --skip_rows 1,2 -i {input[0]} --tree {output[0]} \
            --annotation {output[1]} \
            --most_abundant 100 --abundance_threshold 1 \
            --least_biomarkers 10 --annotations 5,6 \
            --external_annotations 7 --min_clade_size 1
        """

rule graphlan_annotate:
    input:
        opj(config["report_path"],"metaphlan2","metaphlan2.tree.txt"),
        opj(config["report_path"],"metaphlan2","metaphlan2.annot.txt")
    output:
        temp(opj(config["report_path"],"metaphlan2","metaphlan2.abundance.xml"))
    conda: 
        "../../../envs/graphlan.yml"
    shell:
        """
        graphlan_annotate.py --annot {input[1]} {input[0]} {output[0]}
        """

rule graphlan:
    input:
        opj(config["report_path"],"metaphlan2","metaphlan2.abundance.xml")
    output:
        opj(config["report_path"],"metaphlan2","metaphlan2.graph.png")
    conda: 
        "../../../envs/graphlan.yml"
    shell:
        """
        graphlan.py --dpi 300 {input[0]} {output[0]}
        """