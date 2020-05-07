localrules: download_gtdb, aggregate_gtdbtk

rule download_gtdb:
    output:
        met = opj(config["resource_path"], "gtdb", "metadata",
                  "metadata.txt")
    log:
        opj(config["resource_path"], "gtdb", "download.log")
    params:
        url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz",
        tar = lambda w, output: opj(os.path.dirname(output.met), "gtdbtk_r89_data.tar.gz"),
        dir = lambda w, output: os.path.dirname(output.met)
    shell:
        """
        curl -L -o {params.tar} {params.url} > {log} 2>&1
        tar xzf {params.tar} -C {params.dir} --strip 1 > {log} 2>&1
        """

rule gtdbtk_classify:
    input:
        met = opj(config["resource_path"], "gtdb", "metadata",
                  "metadata.txt"),
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv")
    output:
        touch(opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "gtdbtk", "done"))
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "gtdbtk", "gtdbtk.log")
    params:
        suff = 'fa',
        indir = lambda wildcards: get_indir(wildcards),
        dbdir = lambda wildcards, input: os.path.abspath(os.path.dirname(os.path.dirname(input.met))),
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/gtdbtk.yaml"
    shell:
        """
        bins=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            echo "NO BINS FOUND" > {output}
        else
            export PYTHONPATH=$(which python)
            export GTDBTK_DATA_PATH={params.dbdir}
            gtdbtk classify_wf -x {params.suff} --out_dir {params.outdir} \
                --cpus {threads} --pplacer_cpus {threads} \
                --genome_dir {params.indir} > {log} 2>&1
        fi
        """

rule aggregate_gtdbtk:
    """
    Aggregates GTDB-TK phylogenetic results from several assemblies into a 
    single table.
    """
    input:
        expand(opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "gtdbtk", "done"),
               binner = get_binners(config),
               group = assemblyGroups.keys(),
               l = config["min_contig_length"])
    output:
        summary = opj(config["report_path"], "gtdbtk", "gtdbtk.summary.tsv")
    run:
        summaries = []
        for f in input:
            gtdb_dir = os.path.dirname(f)
            for m in ["bac120", "ar122"]:
                summary = opj(gtdb_dir, "gtdbtk.{}.summary.tsv".format(m))
                if os.path.exists(summary):
                    summaries.append(summary)
        df = concatenate(summaries)
        df.to_csv(output.summary, sep="\t")
