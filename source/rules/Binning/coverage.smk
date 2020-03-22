## COVERAGE FOR CONCOCT ##
rule concoct_coverage_table:
    input:
        bam=get_all_files(samples, opj(config["results_path"],"assembly",
                                       "{group}", "mapping"), ".bam"),
        bed=opj(config["results_path"],"assembly","{group}",
                "final_contigs_cutup.bed")
    output:
        cov=opj(config["results_path"],"concoct","{group}",
                "cov","concoct_inputtable.tsv")
    conda:
        "../../../envs/concoct.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    params:
        samplenames=opj(config["results_path"],"concoct",
                        "{group}","cov","samplenames")
    shell:
        """
        for f in {input.bam} ; 
            do 
                n=$(basename $f); 
                s=$(echo -e $n | sed 's/_[ps]e.markdup.bam//g'); 
                echo $s; 
            done > {params.samplenames}
        concoct_coverage_table.py \
            --samplenames {params.samplenames} \
            {input.bed} {input.bam} > {output.cov}
        rm {params.samplenames}
        """

## COVERAGE FOR METABAT2 ##
rule metabat_coverage:
    input:
        bam=get_all_files(samples, opj(config["results_path"],"assembly",
                                       "{group}","mapping"),".bam")
    output:
        depth=opj(config["results_path"],"metabat","{group}","cov","depth.txt")
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    conda:
        "../../../envs/metabat.yml"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} {input.bam} 
        """