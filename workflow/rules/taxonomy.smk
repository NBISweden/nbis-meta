localrules:
    taxonomy,
    assign_orfs,


##### taxonomy master rule #####
rule taxonomy:
    input:
        expand(
            results + "/annotation/{assembly}/{seqTaxDB}.taxonomy.tsv",
            assembly=assemblies.keys(),
            seqTaxDB=config["taxonomy"]["database"],
        ),
        expand(
            results + "/annotation/{assembly}/{seqTaxDB}.orfs.taxonomy.tsv",
            assembly=assemblies.keys(),
            seqTaxDB=config["taxonomy"]["database"],
        ),


rule mmseq_downloadDB:
    output:
        expand(
            "resources/mmseqs/{{seqTaxDB}}{suff}",
            suff=[
                "",
                "_taxonomy",
                "_mapping",
                "_h.dbtype",
                "_h.index",
                "_h",
                ".lookup",
                ".dbtype",
                ".index",
            ],
        ),
    params:
        seqTaxDB="{seqTaxDB}",
        tmpdir="$TMPDIR.{seqTaxDB}",
    threads: 4
    resources:
        runtime=24 * 60 * 3,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "../envs/mmseqs.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        mmseqs databases  {params.seqTaxDB} {output[0]} {params.tmpdir} --threads {threads}
        rm -rf {params.tmpdir}
        """


rule mmseqs_createqueryDB:
    output:
        queryDB=expand(
            results + "/annotation/{{assembly}}/mmseqs/queryDB{suff}",
            suff=["", "_h", ".index", "_h.index"],
        ),
    input:
        fa=results + "/assembly/{assembly}/final_contigs.fa",
    log:
        results + "/annotation/{assembly}/mmseqs/createdb.log",
    params:
        mem_mib=mem_allowed,
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "../envs/mmseqs.yml"
    resources:
        runtime=60 * 2,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    shell:
        """
        mmseqs createdb {input.fa} {output.queryDB[0]} --dbtype 2 --shuffle 0 \
            --createdb-mode 1 --write-lookup 0 --id-offset 0 --compressed 0 -v 3 > {log} 2>&1
        """


rule mmseqs_taxonomy:
    output:
        expand(
            results + "/annotation/{{assembly}}/mmseqs/{{seqTaxDB}}{suff}",
            suff=[".dbtype", ".index", "_aln.dbtype", "_aln.index"],
        ),
    input:
        queryDB=rules.mmseqs_createqueryDB.output.queryDB[0],
        seqTaxDB="resources/mmseqs/{seqTaxDB}",
    log:
        results + "/annotation/{assembly}/mmseqs/{seqTaxDB}.taxonomy.log",
    params:
        lca_ranks=",".join(config["taxonomy"]["ranks"]),
        out=lambda wildcards, output: os.path.dirname(output[0])
        + "/"
        + wildcards.seqTaxDB,
        tmpdir="$TMPDIR/mmseqs.taxonomy.{assembly}",
        mem_mib=mem_allowed,
        extra_params=config["taxonomy"]["mmseqs_extra_params"],
    threads: 20
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "../envs/mmseqs.yml"
    resources:
        runtime=60 * 10,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    shell:
        """
        mkdir -p {params.tmpdir}
        mmseqs taxonomy {input.queryDB} {input.seqTaxDB} {params.out} {params.tmpdir} \
            --lca-mode 3 --tax-output-mode 2 --lca-ranks {params.lca_ranks} --tax-lineage 1 \
            --threads {threads} --local-tmp {params.tmpdir} --remove-tmp-files 1 > {log} 2>&1 
        """


rule mmseqs_createtsv:
    output:
        lca=results + "/annotation/{assembly}/{seqTaxDB}.taxonomy.tsv",
    input:
        queryDB=rules.mmseqs_createqueryDB.output.queryDB[0],
        taxRes=rules.mmseqs_taxonomy.output,
    log:
        results + "/annotation/{assembly}/mmseqs/{seqTaxDB}.createtsv.log",
    params:
        taxRes=lambda wildcards, input: os.path.dirname(input.taxRes[0])
        + "/"
        + wildcards.seqTaxDB,
        mem_mib=mem_allowed,
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "../envs/mmseqs.yml"
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    shell:
        """
        mmseqs createtsv {input.queryDB} {params.taxRes} {output.lca} \
            --threads {threads} --first-seq-as-repr 0 --target-column 1 \
            --full-header 0 --idx-seq-src 0 --db-output 0 --compressed 0 -v 3  
        """


rule assign_orfs:
    input:
        tax=rules.mmseqs_createtsv.output.lca,
        gff=results + "/annotation/{assembly}/final_contigs.gff",
    output:
        tax=results + "/annotation/{assembly}/{seqTaxDB}.orfs.taxonomy.tsv",
    params:
        ranks=config["taxonomy"]["ranks"],
    script:
        "../scripts/taxonomy_utils.py"
