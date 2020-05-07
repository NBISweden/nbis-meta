rule barrnap:
    """
    Identify rRNA genes in genome bins
    """
    input:
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv"),
        gtdbtk = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "gtdbtk", "done")
    output:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "barrnap", "rRNA.gff")
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "barrnap", "log")
    conda:
        "../../../envs/barrnap.yaml"
    params:
        indir = lambda wildcards: get_indir(wildcards),
        gtdbtk_dir = lambda w, input: os.path.dirname(input.gtdbtk)
    resources:
        runtime = lambda wildcards, attempt: attempt**2*30
    threads: 1
    shell:
        """
        bins=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            touch {output}
        else
            cat {params.gtdbtk_dir}/gtdbtk.*.summary.tsv | cut -f1 -d ';' | grep -v "user_genome" | \
                while read line;
                do
                    d=$(echo -e "$line" | cut -f2)
                    g=$(echo -e "$line" | cut -f1)
                    if [ "$d" == "d__Bacteria" ]; then
                        k="bac"
                    else
                        k="arc"
                    fi
                    barrnap --kingdom $k --quiet {params.indir}/$g.fa | \
                        egrep -v "^#" | sed "s/$/;genome=$g/g" >> {output}
                done
        fi              
        """

rule count_rRNA:
    input:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "barrnap", "rRNA.gff")
    output:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "barrnap", "rRNA.types.tsv")
    script:
        "../../../scripts/count_rRNA.py"

rule trnascan:
    """
    Identify tRNA genes in genome bins
    """
    input:
        tsv = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "summary_stats.tsv"),
        gtdbtk = opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "gtdbtk", "done")
    output:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "tRNAscan", "tRNA.tsv")
    log:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "tRNAscan", "tRNA.log")
    params:
        indir = lambda wildcards: get_indir(wildcards),
        gtdbtk_dir = lambda w, input: os.path.dirname(input.gtdbtk)
    resources:
        runtime = lambda wildcards, attempt: attempt*30
    threads: 4
    conda:
        "../../../envs/annotation.yaml"
    shell:
        """
        bins=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $bins == 0 ] ; then
            touch {output}
        else
            echo -e "Name\ttRNA#\ttRNA_Begin\ttRNA_End\ttRNA_type\tAnti_Codon\tIntron_Begin\tIntron_End\tInf_Score\tNote\tBin_Id" > {output}
            cat {params.gtdbtk_dir}/gtdbtk.*.summary.tsv | cut -f1 -d ';' | grep -v "user_genome" | \
                while read line;
                do
                    d=$(echo -e "$line" | cut -f2)
                    g=$(echo -e "$line" | cut -f1)
                    if [ "$d" == "d__Bacteria" ]; then
                        model="-B"
                    else
                        model="-A"
                    fi
                    tRNAscan-SE $model --quiet --thread {threads} \
                        {params.indir}/$g.fa | tail -n +4 | sed "s/$/\t$g/g" >> {output}
                done
        fi          
        """

rule count_tRNA:
    input:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "tRNAscan", "tRNA.tsv")
    output:
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "tRNAscan", "tRNA.types.tsv"),
        opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "tRNAscan", "tRNA.total.tsv")
    script:
        "../../../scripts/count_tRNA.py"

rule aggregate_bin_annot:
    input:
        trna = expand(opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "tRNAscan", "tRNA.total.tsv"),
               binner = get_binners(config),
               group = assemblyGroups.keys(),
               l = config["min_contig_length"]),
        rrna = expand(opj(config["results_path"], "binning", "{binner}", "{group}", "{l}", "barrnap", "rRNA.types.tsv"),
               binner = get_binners(config),
               group = assemblyGroups.keys(),
               l = config["min_contig_length"])
    output
        trna = opj(config["report_path"], "bin_annotation", "tRNA.total.tsv"),
        rrna = opj(config["report_path"], "bin_annotation", "rRNA.types.tsv")
    run:
        df = concatenate(input.trna)
        df.to_csv(output.trna, sep="\t", index=True)
        df = concatenate(input.rrna)
        df.to_csv(output.rrna, sep="\t", index=True)
