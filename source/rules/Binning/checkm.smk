localrules: download_checkm, checkm_qa, aggregate_checkm_stats, checkm_profile

rule download_checkm:
    output:
        db = opj(config["resource_path"], "checkm", ".dmanifest")
    log:
        opj(config["resource_path"], "checkm", "checkm.log")
    params:
        tar = lambda wildcards, output: opj(os.path.dirname(output.db), "checkm_data.tar.gz"),
        dir = lambda wildcards, output: os.path.dirname(output.db)
    conda:
        "../../../envs/checkm.yaml"
    shell:
        """
        # Download
        curl -L https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -o {params.tar} -s
        # Extract
        tar -C {params.dir} -xf {params.tar}
        # Set root
        checkm data setRoot {params.dir} > {log} 2>&1
        """

rule checkm_lineage_wf:
    input:
        db = opj(config["resource_path"], "checkm", ".dmanifest"),
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv")
    output:
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
                  "genome_stats.tsv"),
        ms = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
                  "lineage.ms")
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
            "checkm.log")
    conda:
        "../../../envs/checkm.yaml"
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    params:
        suff = 'fa',
        indir = lambda wildcards: get_indir(wildcards),
        outdir = lambda wildcards, output: os.path.dirname(output.tsv),
        tree = get_tree_settings(config)
    shell:
        """
        bins=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            echo "NO BINS FOUND" > {output.tsv}
            touch {output.ms}
        else
            checkm lineage_wf -t {threads} --pplacer_threads {threads} \
                -x {params.suff} {params.tree} -q \
                --tab_table -f {output.tsv} \
                {params.indir} {params.outdir} \
                > {log} 2>&1
        fi
        """

rule checkm_qa:
    """
    Runs checkm qa to generate output format 2 with extended summaries of bins
    """
    input:
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv"),
        ms = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
                  "lineage.ms")
    output:
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
                  "genome_stats.extended.tsv")
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
                  "qa.log")
    conda:
        "../../../envs/checkm.yaml"
    params:
        dir = lambda wildcards, output: os.path.dirname(output.tsv)
    shell:
        """
        bins=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            echo "NO BINS FOUND" > {output.tsv}
        else
            checkm qa -o 2 --tab_table -f {output.tsv} \
                {input.ms} {params.dir} > {log} 2>&1
        fi 
        """

rule checkm_coverage:
    input:
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv"),
        bam = get_all_files(samples, opj(config["results_path"], "assembly",
                                         "{group}", "mapping"), ".markdup.bam"),
        bai = get_all_files(samples, opj(config["results_path"], "assembly",
                                         "{group}", "mapping"), ".markdup.bam.bai")
    output:
        cov = temp(opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm", "coverage.tsv"))
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm", "checkm_coverage.log")
    params:
        dir = lambda wildcards: get_indir(wildcards)
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/checkm.yaml"
    shell:
        """
        bins=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            echo "NO BINS FOUND" > {output}
        else
            checkm coverage -x fa -t {threads} {params.dir} \
                {output} {input.bam} > {log} 2>&1
        fi
        """

rule checkm_profile:
    input:
        cov = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm", "coverage.tsv"),
        stats = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv")
    output:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm", "profile.tsv")
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm", "checkm_profile.log")
    conda:
        "../../../envs/checkm.yaml"
    shell:
        """
        bins=$(wc -l {input.stats} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            echo "NO BINS FOUND" > {output}
        else
            checkm profile -f {output} --tab_table {input.cov} > {log} 2>&1
        fi
        """

rule aggregate_checkm_profiles:
    input:
        expand(opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "checkm",
                  "profile.tsv"),
               group = assemblyGroups.keys(),
               l = config["min_contig_length"],
               binner = get_binners(config))
    output:
        tsv = opj(config["report_path"], "checkm", "checkm.profiles.tsv")
    run:
        df = concatenate(input)
        df.to_csv(output.tsv, sep="\t", index=True)

rule aggregate_checkm_stats:
    input:
        expand(opj(config["results_path"], "binning", "{binner}", "{group}",
                   "{l}", "checkm", "genome_stats.extended.tsv"),
               group = assemblyGroups.keys(),
               l = config["min_contig_length"],
               binner = get_binners(config))
    output:
        tsv = opj(config["report_path"], "checkm", "checkm.stats.tsv")
    run:
        df = concatenate(input)
        df.to_csv(output.tsv, sep="\t", index=True)
