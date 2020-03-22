localrules:
    concoct_cutup,
    merge_cutup,
    extract_fasta,
    maxbin_stats,
    concoct_stats,
    metabat_stats

## METABAT2 ##

rule run_metabat:
    input:
        fa=opj(config["results_path"],"assembly","{group}",
               "final_contigs.fa"),
        depth=opj(config["results_path"],"metabat","{group}",
                  "cov","depth.txt")
    output:
        opj(config["results_path"],"metabat","{group}","{l}","contig_map.tsv")
    conda:
        "../../../envs/metabat.yml"
    threads: config["metabat_threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    params:
        n=opj(config["results_path"],"metabat","{group}","{l}","metabat")
    shell:
        """
        metabat2 \
            -i {input.fa} \
            -a {input.depth} \
            -m {wildcards.l} \
            -t {threads} \
            -o {params.n}
        grep '>' {params.n}*.fa | \
            awk -F/ '{{print $NF}}' | \
            sed 's/.fa:>/\t/g' > {output[0]}
        """

rule metabat_stats:
    input:
        opj(config["results_path"],"metabat","{group}",
            "{l}","contig_map.tsv")
    output:
        opj(config["results_path"],"metabat","{group}",
            "{l}","summary_stats.tsv")
    params:
        dir=opj(config["results_path"],"metabat","{group}","{l}"),
        suffix=".fa"
    shell:
        """
        python source/utils/binning_stats.py \
            --suffix {params.suffix} {params.dir} > {output[0]}
        """

## MAXBIN ##

def get_fw_reads(config, p):
    """
    MaxBin2 only uses unpaired reads for mapping with bowtie2.
    Here we iterate over all samples
    """
    files=[]
    for sample in samples.keys():
        for run in samples[sample].keys():
            if "R1" in samples[sample][run].keys():
                f=opj(config["intermediate_path"],"preprocess",
                        "{sample}_{run}_R1{p}.fastq.gz".format(sample=sample,
                                                               run=run,
                                                               p=p))
            else:
                f=opj(config["intermediate_path"],"preprocess",
                      "{sample}_{run}_se{p}.fastq.gz".format(sample=sample,
                                                             run=run,
                                                             p=p))
            files.append(f)
    reads_string=""
    for i, f in enumerate(files, start=1):
        reads_string+="-reads{i} {f} ".format(i=i, f=f)
    return reads_string

rule run_maxbin:
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        opj(config["results_path"],"maxbin","{group}","{l}","{group}.summary")
    log:
        opj(config["results_path"],"maxbin","{group}","{l}","maxbin.log")
    params:
        dir=opj(config["results_path"],"maxbin","{group}","{l}"),
        tmp_dir=opj(config["scratch_path"],"{group}","{l}"),
        reads=get_fw_reads(config, PREPROCESS),
        markerset=config["maxbin_markerset"]
    threads: config["maxbin_threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*5
    conda:
        "../../../envs/maxbin.yml"
    shell:
        """
        mkdir -p {params.dir}
        mkdir -p {params.tmp_dir}
        run_MaxBin.pl \
            -markerset {params.markerset} \
            -contig {input} \
            {params.reads} \
            -min_contig_length {wildcards.l} \
            -thread {threads} \
            -out {params.tmp_dir}/{wildcards.group} 2>{log}
        mv {params.tmp_dir}/* {params.dir}
        rm -r {params.tmp_dir}
        """

rule maxbin_stats:
    input:
        opj(config["results_path"],"maxbin","{group}","{l}","{group}.summary")
    output:
        opj(config["results_path"],"maxbin","{group}",
            "{l}","summary_stats.tsv")
    params:
        dir=opj(config["results_path"],"maxbin","{group}","{l}"),
        suffix=".fasta"
    shell:
        """
        python source/utils/binning_stats.py \
            --suffix {params.suffix} {params.dir} > {output[0]}
        """

## CONCOCT ##

rule concoct_cutup:
    input:
        fa=opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        fa=opj(config["results_path"],"assembly","{group}",
               "final_contigs_cutup.fa"),
        bed=opj(config["results_path"],"assembly","{group}",
                "final_contigs_cutup.bed")
    conda:
        "../../../envs/concoct.yml"
    shell:
        """
        cut_up_fasta.py \
            -b {output.bed} -c 10000 -o 0 -m {input.fa} > {output.fa}
        """

rule run_concoct:
    input:
        cov=opj(config["results_path"],"concoct","{group}",
                "cov","concoct_inputtable.tsv"),
        fa=opj(config["results_path"],"assembly","{group}",
               "final_contigs_cutup.fa")
    output:
        opj(config["results_path"],"concoct","{group}","{l}",
            "clustering_gt{l}.csv")
    log:
        opj(config["results_path"],"concoct","{group}","{l}","log.txt")
    params:
        basename=opj(config["results_path"],"concoct","{group}","{l}"),
        length="{l}"
    threads: config["concoct_threads"]
    conda:
        "../../../envs/concoct.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    shell:
        """
        concoct \
            -t {threads} \
            --coverage_file {input.cov} \
            --composition_file {input.fa} \
            -b {params.basename}/ \
            -l {params.length} >/dev/null 2>&1
        """

rule merge_cutup:
    input:
        opj(config["results_path"],"concoct","{group}",
            "{l}","clustering_gt{l}.csv")
    output:
        opj(config["results_path"],"concoct","{group}",
            "{l}","clustering_gt{l}_merged.csv"),
        opj(config["results_path"],"concoct","{group}",
            "{l}","clustering_gt{l}_merged.log")
    conda:
        "../../../envs/concoct.yml"
    shell:
        """
        merge_cutup_clustering.py \
            {input[0]} > {output[0]} 2> {output[1]}
        """

rule extract_fasta:
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa"),
        opj(config["results_path"],"concoct","{group}",
            "{l}","clustering_gt{l}_merged.csv")
    output:
        touch(opj(config["results_path"],"concoct","{group}",
                  "{l}","fasta","done"))
    params:
        dir=opj(config["results_path"],"concoct","{group}","{l}","fasta")
    conda:
        "../../../envs/concoct.yml"
    shell:
        """
        extract_fasta_bins.py {input[0]} {input[1]} --output_path {params.dir}
        """

rule concoct_stats:
    input:
        opj(config["results_path"],"concoct","{group}","{l}","fasta","done")
    output:
        opj(config["results_path"],"concoct","{group}",
            "{l}","summary_stats.tsv")
    params:
        dir=opj(config["results_path"],"concoct","{group}","{l}","fasta"),
        suffix=".fa"
    shell:
        """
        python source/utils/binning_stats.py \
            --suffix {params.suffix} {params.dir} > {output[0]}
        """