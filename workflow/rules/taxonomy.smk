localrules:
    taxonomy,
    download_sourmash_db,
    contigtax_assign_orfs,
    sourmash_compute,
    merge_contigtax_sourmash,


##### taxonomy master rule #####
rule taxonomy:
    input:
        expand(
            results + "/annotation/{assembly}/taxonomy/orfs.{db}.taxonomy.tsv",
            assembly=assemblies.keys(),
            db=config["taxonomy"]["database"],
        ),


rule mmseqs_taxonomy:
    output:
        tsv = results + "/annotation/{assembly}/{seqTaxDB}_lca.tsv",
    log:
        results + "/annotation/{assembly}/{seqTaxDB}.mmseqs.log"
    input:
        seqTaxDB=expand(
            "resources/mmseqs/{{seqTaxDB}}{suff}",
            suff=[
                "",
                ".index",
                ".dbtype",
                ".lookup",
                "_h",
                "_h.dbtype",
                "_h.index",
                "_mapping",
                "_taxonomy",
            ],
        ),
        fa = results + "/assembly/{assembly}/final_contigs.fa"
    params:
        seqTaxDB="resources/mmseqs/{seqTaxDB}",
        tmpdir="$TMPDIR/mmseqs.{assembly}",
        out = lambda wildcards, output: os.path.dirname(output.tsv) + "/" + wildcards.seqTaxDB,
        lca_ranks = ",".join(config["taxonomy"]["ranks"]),
        sensitivity = 7.5
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "../envs/mmseqs.yml"
    threads: 10
    shell:
        """
        rm -rf {params.tmpdir}*
        mkdir -p {params.tmpdir}
        mmseqs easy-taxonomy --local-tmp {params.tmpdir} --lca-mode 3 --tax-lineage 1 \
            --lca-ranks {params.lca_ranks} --threads {threads} \
            {input.fa} {params.seqTaxDB} {params.out} {params.tmpdir} > {log} 2>&1
        rm -rf {params.tmpdir}*
        """

##### sourmash #####


rule download_sourmash_db:
    output:
        lca="resources/sourmash/sourmash_db.lca.json",
        version="resources/sourmash/version.txt",
    log:
        "resources/sourmash/download.log",
    params:
        url=config["taxonomy"]["sourmash_database_url"],
    shell:
        """
        curl -L -v -o {output.lca}.gz {params.url} > {log} 2>&1
        grep filename {log} | cut -f2 -d ';' > {output.version}
        gunzip {output.lca}.gz
        """


rule sourmash_compute:
    input:
        results + "/assembly/{assembly}/final_contigs.fa",
    output:
        results + "/assembly/{assembly}/final_contigs.fa.sig",
    log:
        results + "/assembly/{assembly}/sourmash_compute.log",
    conda:
        "../envs/sourmash.yml"
    params:
        frac=config["taxonomy"]["sourmash_fraction"],
        k=config["taxonomy"]["sourmash_kmer_size"],
    shell:
        """
        sourmash compute --singleton --scaled {params.frac} \
            -k {params.k} -o {output} {input} > {log} 2>&1
        """


rule sourmash_classify:
    input:
        sig=results + "/assembly/{assembly}/final_contigs.fa.sig",
        db="resources/sourmash/sourmash_db.lca.json",
    output:
        csv=results + "/annotation/{assembly}/taxonomy/sourmash.taxonomy.csv",
    log:
        results + "/annotation/{assembly}/taxonomy/sourmash.log",
    params:
        frac=config["taxonomy"]["sourmash_fraction"],
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2 * 30,
    conda:
        "../envs/sourmash.yml"
    shell:
        """
        sourmash lca classify --db {input.db} --scaled {params.frac} \
            --query {input.sig} -o {output.csv} > {log} 2>&1
        """


##### common taxonomy rules #####


rule merge_contigtax_sourmash:
    input:
        smash=results + "/annotation/{assembly}/taxonomy/sourmash.taxonomy.csv",
        contigtax=expand(
            results + "/annotation/{{assembly}}/taxonomy/contigtax.{db}.taxonomy.tsv",
            db=config["taxonomy"]["database"],
        ),
    output:
        results + "/annotation/{assembly}/taxonomy/final_contigs.taxonomy.tsv",
    log:
        results + "/annotation/{assembly}/taxonomy/merge.log",
    script:
        "../scripts/taxonomy_utils.py"


rule contigtax_assign_orfs:
    input:
        tax=results + "/annotation/{assembly}/taxonomy/final_contigs.taxonomy.tsv",
        gff=results + "/annotation/{assembly}/final_contigs.gff",
    output:
        tax=expand(
            results + "/annotation/{{assembly}}/taxonomy/orfs.{db}.taxonomy.tsv",
            db=config["taxonomy"]["database"],
        ),
    script:
        "../scripts/taxonomy_utils.py"
