from scripts.common import binning_input, get_fw_reads
from scripts.common import get_binners, get_tree_settings, concatenate

localrules:
    bin,
    concoct_cutup,
    merge_cutup,
    extract_fasta,
    contig_map,
    binning_stats,
    aggregate_binning_stats,
    download_checkm,
    checkm_qa,
    remove_checkm_zerocols,
    aggregate_checkm_stats,
    aggregate_checkm_profiles,
    checkm_profile,
    download_gtdb,
    aggregate_gtdbtk,
    download_ref_genome,
    generate_fastANI_lists,
    cluster_genomes,
    count_rRNA,
    count_tRNA,
    aggregate_bin_annot,
    binning_report

##### master rule for binning #####

rule bin:
    input:
        binning_input(config),
        results+"/report/binning/bin_report.pdf"

##### target rule for running checkm analysis #####

rule checkm:
    input:
        results+"/report/checkm/checkm.stats.tsv"

##### metabat2 #####

rule metabat_coverage:
    input:
        bam=get_all_files(samples, results+"/assembly/{assembly}/mapping", ".bam")
    output:
        depth=results+"/binning/metabat/{assembly}/cov/depth.txt"
    log:
        results+"/binning/metabat/{assembly}/cov/log"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    conda:
        "../envs/metabat.yml"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} {input.bam} >{log} 2>&1
        """

rule metabat:
    input:
        fa=results+"/assembly/{assembly}/final_contigs.fa",
        depth=results+"/binning/metabat/{assembly}/cov/depth.txt"
    output:
        touch(results+"/binning/metabat/{assembly}/{l}/done")
    log:
        results+"/binning/metabat/{assembly}/{l}/metabat.log"
    params:
        n=results+"/binning/metabat/{assembly}/{l}/metabat"
    conda:
        "../envs/metabat.yml"
    threads: config["binning"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        metabat2 -i {input.fa} -a {input.depth} -m {wildcards.l} -t {threads} \
            -o {params.n} > {log} 2>&1
        """

##### maxbin2 #####

rule maxbin:
    input:
        results+"/assembly/{assembly}/final_contigs.fa"
    output:
        touch(results+"/binning/maxbin/{assembly}/{l}/done")
    log:
        results+"/binning/maxbin/{assembly}/{l}/maxbin.log"
    params:
        dir=results+"/binning/maxbin/{assembly}/{l}",
        tmp_dir=temppath+"/maxbin/{assembly}/{l}",
        reads=get_fw_reads(config, samples, PREPROCESS),
        markerset=config["maxbin"]["markerset"]
    threads: config["binning"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*5
    conda:
        "../envs/maxbin.yml"
    shell:
        """
        set +e
        mkdir -p {params.dir}
        mkdir -p {params.tmp_dir}
        run_MaxBin.pl -markerset {params.markerset} -contig {input} \
            {params.reads} -min_contig_length {wildcards.l} -thread {threads} \
            -out {params.tmp_dir}/maxbin >{log} 2>{log}
        exitcode=$?
        if [ $exitcode -eq 255 ]; then
            exit 0
        else
            # Rename fasta files
            ls {params.tmp_dir} | grep ".fasta" | while read f;
            do
                mv {params.tmp_dir}/$f {params.dir}/${{f%.fasta}}.fa
            done
            # Move output from temporary dir
            ls {params.tmp_dir} | while read f;
            do
                mv {params.tmp_dir}/$f {params.dir}/
            done
        fi
        # Clean up
        rm -r {params.tmp_dir}
        """

##### concoct #####

rule concoct_coverage_table:
    input:
        bam=get_all_files(samples, results+"/assembly/{assembly}/mapping", ".bam"),
        bai=get_all_files(samples, results+"/assembly/{assembly}/mapping", ".bam.bai"),
        bed=results+"/assembly/{assembly}/final_contigs_cutup.bed"
    output:
        cov=results+"/binning/concoct/{assembly}/cov/concoct_inputtable.tsv"
    conda:
        "../envs/concoct.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    params:
        samplenames=results+"/binning/concoct/{assembly}/cov/samplenames",
        p=POSTPROCESS
    shell:
        """
        for f in {input.bam} ;
            do
                n=$(basename $f);
                s=$(echo -e $n | sed 's/_[ps]e{params.p}.bam//g');
                echo $s;
            done > {params.samplenames}
        concoct_coverage_table.py \
            --samplenames {params.samplenames} \
            {input.bed} {input.bam} > {output.cov}
        rm {params.samplenames}
        """

rule concoct_cutup:
    input:
        fa=results+"/assembly/{assembly}/final_contigs.fa"
    output:
        fa=results+"/assembly/{assembly}/final_contigs_cutup.fa",
        bed=results+"/assembly/{assembly}/final_contigs_cutup.bed"
    log:
        results+"/assembly/{assembly}/final_contigs_cutup.log"
    conda:
        "../envs/concoct.yml"
    shell:
        """
        cut_up_fasta.py -b {output.bed} -c 10000 -o 0 -m {input.fa} \
            > {output.fa} 2>{log}
        """

rule concoct:
    input:
        cov=results+"/binning/concoct/{assembly}/cov/concoct_inputtable.tsv",
        fa=results+"/assembly/{assembly}/final_contigs_cutup.fa"
    output:
        results+"/binning/concoct/{assembly}/{l}/clustering_gt{l}.csv"
    log:
        results+"/binning/concoct/{assembly}/{l}/log.txt"
    params:
        basename=lambda wildcards, output: os.path.dirname(output[0]),
        length="{l}"
    threads: config["binning"]["threads"]
    conda:
        "../envs/concoct.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    shell:
        """
        concoct -t {threads} --coverage_file {input.cov} -l {params.length} \
            --composition_file {input.fa} -b {params.basename}/ >/dev/null 2>&1
        """

rule merge_cutup:
    input:
        results+"/binning/concoct/{assembly}/{l}/clustering_gt{l}.csv"
    output:
        results+"/binning/concoct/{assembly}/{l}/clustering_gt{l}_merged.csv"
    log:
        results+"/binning/concoct/{assembly}/{l}/clustering_gt{l}_merged.log"
    conda:
        "../envs/concoct.yml"
    shell:
        """
        merge_cutup_clustering.py {input[0]} > {output[0]} 2> {log}
        """

rule extract_fasta:
    input:
        results+"/assembly/{assembly}/final_contigs.fa",
        results+"/binning/concoct/{assembly}/{l}/clustering_gt{l}_merged.csv"
    output:
        touch(results+"/binning/concoct/{assembly}/{l}/done")
    log:
        results+"/binning/concoct/{assembly}/{l}/extract_fasta.log"
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0]),
        tmp_dir=temppath+"/concoct/{assembly}/{l}"
    conda:
        "../envs/concoct.yml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        extract_fasta_bins.py {input[0]} {input[1]} \
            --output_path {params.tmp_dir} 2> {log}
        ls {params.tmp_dir} | egrep "[0-9].fa" | while read f;
        do
            mv {params.tmp_dir}/$f {params.dir}/concoct.$f
        done
        """

##### map contigs to bins #####

rule contig_map:
    input:
        results+"/binning/{binner}/{assembly}/{l}/done"
    output:
        results+"/binning/{binner}/{assembly}/{l}/contig_map.tsv"
    params:
        dir=lambda wildcards, input: os.path.dirname(input[0])
    script:
        "../scripts/binning_utils.py"

##### bin qc #####

rule binning_stats:
    input:
        results+"/binning/{binner}/{assembly}/{l}/contig_map.tsv"
    output:
        results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv"
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0])
    script:
        "../scripts/binning_utils.py"

rule aggregate_binning_stats:
    input:
        expand(results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv",
               assembly=assemblies.keys(),
               l=config["binning"]["contig_lengths"],
               binner=get_binners(config))
    message:
        "Aggregating statistics on binned genomes"
    output:
        report(results+"/report/binning/binning_summary.tsv",
               category="Binning", caption="../report/binning.rst")
    run:
        df=concatenate(input, index=-2)
        df.to_csv(output[0], sep="\t", index=True)

##### checkm #####

rule download_checkm:
    output:
        db="resources/checkm/.dmanifest"
    log:
        "resources/checkm/checkm.log"
    params:
        tar=lambda wildcards, output: os.path.dirname(output.db)+"/checkm_data.tar.gz",
        dir=lambda wildcards, output: os.path.dirname(output.db)
    conda:
        "../envs/checkm.yml"
    shell:
        """
        # Download
        curl -L https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -o {params.tar} -s
        # Extract
        tar -C {params.dir} -xf {params.tar}
        # Set root
        checkm data setRoot {params.dir} > {log} 2>&1
        """

if config["checkm"]["taxonomy_wf"]:
    rule checkm_taxonomy_wf:
        input:
            db="resources/checkm/.dmanifest",
            tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv"
        output:
            tsv=results+"/binning/{binner}/{assembly}/{l}/checkm/genome_stats.tsv",
            ms=results+"/binning/{binner}/{assembly}/{l}/checkm/lineage.ms"
        log:
            results+"/binning/{binner}/{assembly}/{l}/checkm/checkm.log"
        conda:
            "../envs/checkm.yml"
        threads: 10
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60
        params:
            suff='fa',
            indir=lambda wildcards, input: os.path.dirname(input.tsv),
            outdir=lambda wildcards, output: os.path.dirname(output.tsv),
            rank=config["checkm"]["rank"],
            taxon=config["checkm"]["taxon"]
        shell:
            """
            lines=$(wc -l {input.tsv} | cut -f1 -d ' ')
            if [ $lines == 1 ] ; then
                echo "NO BINS FOUND" > {output.tsv}
                touch {output.ms}
            else
                checkm taxonomy_wf -t {threads} -x {params.suff} -q \
                    --tab_table -f {output.tsv} \
                    {params.rank} {params.taxon} {params.indir} {params.outdir} \
                    > {log} 2>&1
                ln -s {params.taxon}.ms {output.ms}
            fi
            """
else:
    rule checkm_lineage_wf:
        input:
            db="resources/checkm/.dmanifest",
            tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv"
        output:
            tsv=results+"/binning/{binner}/{assembly}/{l}/checkm/genome_stats.tsv",
            ms=results+"/binning/{binner}/{assembly}/{l}/checkm/lineage.ms"
        log:
            results+"/binning/{binner}/{assembly}/{l}/checkm/checkm.log"
        conda:
            "../envs/checkm.yml"
        threads: 10
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60
        params:
            suff='fa',
            indir=lambda wildcards, input: os.path.dirname(input.tsv),
            outdir=lambda wildcards, output: os.path.dirname(output.tsv),
            tree=get_tree_settings(config)
        shell:
            """
            lines=$(wc -l {input.tsv} | cut -f1 -d ' ')
            if [ $lines == 0 ] ; then
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
        tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv",
        ms=results+"/binning/{binner}/{assembly}/{l}/checkm/lineage.ms"
    output:
        tsv=results+"/binning/{binner}/{assembly}/{l}/checkm/genome_stats.extended.tsv"
    log:
        results+"/binning/{binner}/{assembly}/{l}/checkm/qa.log"
    conda:
        "../envs/checkm.yml"
    params:
        dir=lambda wildcards, output: os.path.dirname(output.tsv)
    shell:
        """
        lines=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $lines == 1 ] ; then
            echo "NO BINS FOUND" > {output.tsv}
        else
            checkm qa -o 2 --tab_table -f {output.tsv} \
                {input.ms} {params.dir} > {log} 2>&1
        fi
        """

rule checkm_coverage:
    input:
        tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv",
        bam=get_all_files(samples, results+"/assembly/{assembly}/mapping", "{p}.bam".format(p=POSTPROCESS)),
        bai=get_all_files(samples, results+"/assembly/{assembly}/mapping", "{p}.bam.bai".format(p=POSTPROCESS))
    output:
        cov=temp(results+"/binning/{binner}/{assembly}/{l}/checkm/coverage.tsv")
    log:
        results+"/binning/{binner}/{assembly}/{l}/checkm/checkm_coverage.log"
    params:
        dir=lambda wildcards, input: os.path.dirname(input.tsv)
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/checkm.yml"
    shell:
        """
        lines=$(wc -l {input.tsv} | cut -f1 -d ' ')
        if [ $lines == 1 ] ; then
            echo "NO BINS FOUND" > {output}
        else
            checkm coverage -x fa -t {threads} {params.dir} \
                {output} {input.bam} > {log} 2>&1
        fi
        """

rule remove_checkm_zerocols:
    """
    Pre-checks the checkm coverage file and removes zero count bam columns as
    these can generate ZeroDivisionError in downstream rules.
    """
    input:
        cov=results+"/binning/{binner}/{assembly}/{l}/checkm/coverage.tsv"
    output:
        cov=temp(results+"/binning/{binner}/{assembly}/{l}/checkm/_coverage.tsv")
    script:
        "../scripts/binning_utils.py"

rule checkm_profile:
    input:
        cov=results+"/binning/{binner}/{assembly}/{l}/checkm/_coverage.tsv",
        stats=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv"
    output:
        results+"/binning/{binner}/{assembly}/{l}/checkm/profile.tsv"
    log:
        results+"/binning/{binner}/{assembly}/{l}/checkm/checkm_profile.log"
    conda:
        "../envs/checkm.yml"
    shell:
        """
        lines=$(wc -l {input.stats} | cut -f1 -d ' ')
        cov_lines=$(wc -l {input.cov} | cut -f1 -d ' ')
        if [ $lines == 1 ] ; then
            echo "NO BINS FOUND" > {output}
        elif [ $cov_lines == 0 ] ; then
            echo "NO READS MAPPED" > {output}
        else
            checkm profile -f {output} --tab_table {input.cov} > {log} 2>&1
        fi
        """

rule aggregate_checkm_profiles:
    input:
        expand(results+"/binning/{binner}/{assembly}/{l}/checkm/profile.tsv",
               assembly=assemblies.keys(),
               l=config["binning"]["contig_lengths"],
               binner=get_binners(config))
    output:
        tsv=results+"/report/checkm/checkm.profiles.tsv"
    run:
        df=concatenate(input, index=-3)
        df.to_csv(output.tsv, sep="\t", index=True)

rule aggregate_checkm_stats:
    input:
        expand(results+"/binning/{binner}/{assembly}/{l}/checkm/genome_stats.extended.tsv",
               assembly=assemblies.keys(),
               l=config["binning"]["contig_lengths"],
               binner=get_binners(config))
    output:
        tsv=results+"/report/checkm/checkm.stats.tsv"
    run:
        df=concatenate(input, index=-3)
        df.to_csv(output.tsv, sep="\t", index=True)

##### classify bins with gtdb-tk #####

rule download_gtdb:
    output:
        met="resources/gtdb/metadata/metadata.txt"
    log:
        "resources/gtdb/download.log"
    params:
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz",
        tar=lambda wildcards, output: os.path.dirname(output.met)+"/gtdbtk_r89_data.tar.gz",
        dir=lambda wildcards, output: os.path.dirname(output.met)
    shell:
        """
        curl -L -o {params.tar} {params.url} > {log} 2>&1
        tar xzf {params.tar} -C {params.dir} --strip 1 > {log} 2>&1
        """

rule gtdbtk_classify:
    input:
        met="resources/gtdb/metadata/metadata.txt",
        tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv"
    output:
        touch(results+"/binning/{binner}/{assembly}/{l}/gtdbtk/done")
    log:
        results+"/binning/{binner}/{assembly}/{l}/gtdbtk/gtdbtk.log"
    params:
        suff='fa',
        indir=lambda wildcards, input: os.path.dirname(input.tsv),
        dbdir=lambda wildcards, input: os.path.abspath(os.path.dirname(os.path.dirname(input.met))),
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    conda:
        "../envs/gtdbtk.yml"
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
        expand(results+"/binning/{binner}/{assembly}/{l}/gtdbtk/done",
               binner=get_binners(config),
               assembly=assemblies.keys(),
               l=config["binning"]["contig_lengths"])
    output:
        summary=results+"/report/gtdbtk/gtdbtk.summary.tsv"
    run:
        summaries=[]
        for f in input:
            gtdb_dir=os.path.dirname(f)
            for m in ["bac120/ar122"]:
                summary=gtdb_dir+"/gtdbtk.{}.summary.tsv".format(m)
                if os.path.exists(summary):
                    summaries.append(summary)
        df=concatenate(summaries, index=-3)
        df.to_csv(output.summary, sep="\t")

##### annotate bins #####

rule barrnap:
    """
    Identify rRNA genes in genome bins
    """
    input:
        tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv",
        gtdbtk=results+"/binning/{binner}/{assembly}/{l}/gtdbtk/done"
    output:
        results+"/binning/{binner}/{assembly}/{l}/barrnap/rRNA.gff"
    log:
        results+"/binning/{binner}/{assembly}/{l}/barrnap/log"
    conda:
        "../envs/barrnap.yml"
    params:
        indir=lambda wildcards, input: os.path.dirname(input.tsv),
        gtdbtk_dir=lambda wildcards, input: os.path.dirname(input.gtdbtk),
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
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
                    barrnap --kingdom $k {params.indir}/$g.fa > {params.outdir}/out 2>>{log}
                    lines=$(wc -l {params.outdir}/out | cut -f1 -d ' ')
                    if [ $lines -gt 1 ]; then
                        egrep -v "^#" {params.outdir}/out | sed "s/$/;genome=$g/g" >> {output}
                    else
                        touch {output}
                    fi
                done
            rm {params.outdir}/out
        fi
        """

rule count_rRNA:
    input:
        results+"/binning/{binner}/{assembly}/{l}/barrnap/rRNA.gff"
    output:
        results+"/binning/{binner}/{assembly}/{l}/barrnap/rRNA.types.tsv"
    script:
        "../scripts/binning_utils.py"

rule trnascan_bins:
    #TODO: Run with general model if neither bacteria nor archaea
    """
    Identify tRNA genes in genome bins
    """
    input:
        tsv=results+"/binning/{binner}/{assembly}/{l}/summary_stats.tsv",
        gtdbtk=results+"/binning/{binner}/{assembly}/{l}/gtdbtk/done"
    output:
        results+"/binning/{binner}/{assembly}/{l}/tRNAscan/tRNA.tsv"
    log:
        results+"/binning/{binner}/{assembly}/{l}/tRNAscan/tRNA.log"
    params:
        indir=lambda wildcards, input: os.path.dirname(input.tsv),
        gtdbtk_dir=lambda wildcards, input: os.path.dirname(input.gtdbtk)
    resources:
        runtime=lambda wildcards, attempt: attempt*30
    threads: 4
    conda:
        "../envs/annotation.yml"
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
        results+"/binning/{binner}/{assembly}/{l}/tRNAscan/tRNA.tsv"
    output:
        results+"/binning/{binner}/{assembly}/{l}/tRNAscan/tRNA.types.tsv",
        results+"/binning/{binner}/{assembly}/{l}/tRNAscan/tRNA.total.tsv"
    script:
        "../scripts/binning_utils.py"

rule aggregate_bin_annot:
    input:
        trna=expand(results+"/binning/{binner}/{assembly}/{l}/tRNAscan/tRNA.total.tsv",
                    binner=get_binners(config),
                    assembly=assemblies.keys(),
                    l=config["binning"]["contig_lengths"]),
        rrna=expand(results+"/binning/{binner}/{assembly}/{l}/barrnap/rRNA.types.tsv",
                    binner=get_binners(config),
                    assembly=assemblies.keys(),
                    l=config["binning"]["contig_lengths"])
    output:
        trna=results+"/report/bin_annotation/tRNA.total.tsv",
        rrna=results+"/report/bin_annotation/rRNA.types.tsv"
    run:
        df=concatenate(input.trna, index=-3)
        df.to_csv(output.trna, sep="\t", index=True)
        df=concatenate(input.rrna, index=-3)
        df.to_csv(output.rrna, sep="\t", index=True)

##### genome clustering #####

rule download_ref_genome:
    output:
        "resources/ref_genomes/{genome_id}.fna"
    params:
        ftp_base = lambda wildcards: config["fastani"]["ref_genomes"][wildcards.genome_id]
    script:
        "../scripts/binning_utils.py"

rule generate_fastANI_lists:
    input:
        bins=expand(results+"/binning/{binner}/{assembly}/{l}/checkm/genome_stats.extended.tsv",
                    binner = get_binners(config), assembly = assemblies.keys(),
                    l = config["binning"]["contig_lengths"]),
        refs=expand("resources/ref_genomes/{genome_id}.fna",
                    genome_id = config["fastani"]["ref_genomes"].keys())
    output:
        temp(results+"/binning/fastANI/refList"),
        temp(results+"/binning/fastANI/queryList")
    params:
        outdir = lambda wildcards, output: os.path.abspath(os.path.dirname(output[0])),
        completeness = config["fastani"]["min_completeness"],
        contamination = config["fastani"]["max_contamination"]
    message:
        "Generating input lists for fastANI"
    script:
        "../scripts/binning_utils.py"

rule fastANI:
    input:
        results+"/binning/fastANI/refList",
        results+"/binning/fastANI/queryList"
    output:
        results+"/binning/fastANI/out.txt",
        results+"/binning/fastANI/out.txt.matrix"
    log:
        results+"/binning/fastANI/log"
    threads: 8
    params:
        k = config["fastani"]["kmer_size"],
        frag_len = config["fastani"]["frag_len"],
        fraction = config["fastani"]["fraction"],
        indir = lambda wildcards, input: os.path.dirname(input[0])
    conda:
        "../envs/fastani.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    shell:
        """
        fastANI --rl {input[0]} --ql {input[1]} -k {params.k} -t {threads} \
            --fragLen {params.frag_len} --minFraction {params.fraction} \
            --matrix -o {output[0]} > {log} 2>&1
        """

rule cluster_genomes:
    input:
        mat=results+"/binning/fastANI/out.txt.matrix",
        txt=results+"/binning/fastANI/out.txt"
    output:
        results+"/report/binning/genome_clusters.tsv"
    conda:
        "../envs/fastani.yml"
    params:
        thresh = config["fastani"]["threshold"],
        minfrags = config["fastani"]["minfrags"]
    script:
        "../scripts/binning_utils.py"

##### rule to generate summary plots

rule binning_report:
    input:
        binning_input(config, report=True)
    output:
        report(results+"/report/binning/bin_report.pdf",
               category="Binning", caption="../report/binning.rst")
    message:
        "Plot summary stats of binned genomes"
    conda:
        "../envs/plotting.yml"
    notebook:
        "../notebooks/binning_report.py.ipynb"
