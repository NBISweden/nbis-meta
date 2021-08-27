localrules:
    taxonomy,
    contigtax_download,
    contigtax_download_taxonomy,
    contigtax_download_nr_idmap,
    contigtax_format_nr,
    contigtax_format_uniref,
    contigtax_update,
    download_sourmash_db,
    contigtax_assign_orfs,
    sourmash_compute,
    merge_contigtax_sourmash

##### taxonomy master rule #####
rule taxonomy:
    input:
        expand(results+"/annotation/{assembly}/taxonomy/orfs.{db}.taxonomy.tsv",
               assembly=assemblies.keys(), db=config["taxonomy"]["database"])

##### contigtax #####

rule contigtax_download_taxonomy:
    output:
        sqlite="resources/taxonomy/taxonomy.sqlite",
        taxdump="resources/taxonomy/taxdump.tar.gz",
        nodes="resources/taxonomy/nodes.dmp",
        names="resources/taxonomy/names.dmp",
        pkl="resources/taxonomy/taxonomy.sqlite.traverse.pkl"
    log:
        "resources/taxonomy/contigtax.log"
    params:
        taxdir=lambda wildcards, output: os.path.dirname(output.sqlite)
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax download taxonomy -t {params.taxdir} >{log} 2>&1
        """

rule contigtax_download:
    output:
        fasta=temp("resources/{db}/{db}.fasta.gz")
    log:
        "resources/{db}/contigtax_download.log"
    params:
        dldir=lambda wildcards, output: os.path.dirname(output.fasta),
        tmpdir="$TMPDIR"
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax download {wildcards.db} --tmpdir {params.tmpdir} \
            -d {params.dldir} --skip_idmap >{log} 2>{log}
        """

rule contigtax_download_nr_idmap:
    output:
        idmap="resources/nr/prot.accession2taxid.gz"
    log:
        "resources/nr/contigtax_download_idmap.log"
    params:
        dldir=lambda wildcards, output: os.path.dirname(output.idmap)
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax download idmap -d {params.dldir} > {log} 2>&1
        """

rule contigtax_format_uniref:
    input:
        fasta="resources/{db}/{db}.fasta.gz"
    output:
        fasta="resources/{db}/{db}.reformat.fasta.gz",
        idmap="resources/{db}/prot.accession2taxid.gz"
    log:
        "resources/{db}/contigtax_format.log"
    params:
        tmpdir=temppath
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax format -m {output.idmap} --tmpdir {params.tmpdir} \
            {input.fasta} {output.fasta} > {log} 2>&1
        """

rule contigtax_format_nr:
    input:
        fasta="resources/nr/nr.fasta.gz"
    output:
        fasta="resources/nr/nr.reformat.fasta.gz"
    log:
        "resources/nr/contigtax_format.log"
    params:
        tmpdir=temppath
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax format --tmpdir {params.tmpdir} {input.fasta} \
            {output.fasta} > {log} 2>&1
        """

rule contigtax_update:
    input:
        idmap="resources/{db}/prot.accession2taxid.gz"
    output:
        idmap="resources/{db}/prot.accession2taxid.update.gz"
    log:
        "resources/{db}/contigtax_update.log"
    conda:
        "../envs/taxonomy.yml"
    params:
        dir=lambda wildcards, output: os.path.dirname(output.idmap)
    shell:
        """
        # If an idmap file is available, use it to create an updated idmap file
        if [ -e {params.dir}/idmap.tsv.gz ] ; then
            contigtax update {input.idmap} {params.dir}/idmap.tsg.gz \
                {output.idmap} > {log} 2>&1
        # Otherwise, just create a symlink
        else
            cd {params.dir}
            ln -s $(basename {input.idmap}) $(basename {output.idmap})
        fi
        """

rule contigtax_build:
    input:
        fasta="resources/{db}/{db}.reformat.fasta.gz",
        nodes="resources/taxonomy/nodes.dmp",
        idmap="resources/{db}/prot.accession2taxid.update.gz"
    output:
        "resources/{db}/diamond.dmnd"
    log:
        "resources/{db}/diamond.log"
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/taxonomy.yml"
    shell:
         """
         contigtax build -d {output} -p {threads} {input.fasta} \
            {input.idmap} {input.nodes} >{log} 2>&1
         """

rule contigtax_search:
    input:
        db=expand("resources/{db}/diamond.dmnd", db=config["taxonomy"]["database"]),
        fasta=results+"/assembly/{assembly}/final_contigs.fa"
    output:
        expand(results+"/annotation/{{assembly}}/final_contigs.{db}.tsv.gz",
            db=config["taxonomy"]["database"])
    log:
        results+"/annotation/{assembly}/contigtax_search.log"
    params:
        tmpdir=temppath,
        min_len=config["taxonomy"]["min_len"],
        settings=config["taxonomy"]["search_params"]
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/taxonomy.yml"
    shell:
        """
        contigtax search {params.settings} -p {threads} \
            --tmpdir {params.tmpdir} -l {params.min_len} \
            {input.fasta} {input.db} {output} >{log} 2>&1
        """

rule contigtax_assign:
    input:
        tsv=expand(results+"/annotation/{{assembly}}/final_contigs.{db}.tsv.gz",
                    db=config["taxonomy"]["database"]),
        sql=ancient("resources/taxonomy/taxonomy.sqlite")
    output:
        expand(results+"/annotation/{{assembly}}/taxonomy/contigtax.{db}.taxonomy.tsv",
            db=config["taxonomy"]["database"])
    log:
        results+"/annotation/{assembly}/taxonomy/contigtax_assign.log"
    params:
        taxonomy_ranks=" ".join(config["taxonomy"]["ranks"]),
        taxdir="resources/taxonomy",
        settings=config["taxonomy"]["assign_params"]
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../envs/taxonomy.yml"
    shell:
         """
         contigtax assign {params.settings} -p {threads} -m rank_lca \
            --reportranks {params.taxonomy_ranks} -t {params.taxdir} \
            {input.tsv} {output} > {log} 2>&1
         """

##### sourmash #####

rule download_sourmash_db:
    output:
        lca="resources/sourmash/sourmash_db.lca.json",
        version="resources/sourmash/version.txt"
    log:
        "resources/sourmash/download.log"
    params:
        url=config["taxonomy"]["sourmash_database_url"]
    shell:
        """
        curl -L -v -o {output.lca}.gz {params.url} > {log} 2>&1
        grep filename {log} | cut -f2 -d ';' > {output.version}
        gunzip {output.lca}.gz
        """

rule sourmash_compute:
    input:
        results+"/assembly/{assembly}/final_contigs.fa"
    output:
        results+"/assembly/{assembly}/final_contigs.fa.sig"
    log:
        results+"/assembly/{assembly}/sourmash_compute.log"
    conda:
        "../envs/sourmash.yml"
    params:
        frac=config["taxonomy"]["sourmash_fraction"],
        k=config["taxonomy"]["sourmash_kmer_size"]
    shell:
        """
        sourmash compute --singleton --scaled {params.frac} \
            -k {params.k} -o {output} {input} > {log} 2>&1
        """

rule sourmash_classify:
    input:
        sig=results+"/assembly/{assembly}/final_contigs.fa.sig",
        db="resources/sourmash/sourmash_db.lca.json"
    output:
        csv=results+"/annotation/{assembly}/taxonomy/sourmash.taxonomy.csv"
    log:
        results+"/annotation/{assembly}/taxonomy/sourmash.log"
    params:
        frac=config["taxonomy"]["sourmash_fraction"]
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
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
        smash=results+"/annotation/{assembly}/taxonomy/sourmash.taxonomy.csv",
        contigtax=expand(results+"/annotation/{{assembly}}/taxonomy/contigtax.{db}.taxonomy.tsv",
                        db=config["taxonomy"]["database"])
    output:
        results+"/annotation/{assembly}/taxonomy/final_contigs.taxonomy.tsv"
    log:
        results+"/annotation/{assembly}/taxonomy/merge.log"
    script:
        "../scripts/taxonomy_utils.py"

rule contigtax_assign_orfs:
    input:
        tax=results+"/annotation/{assembly}/taxonomy/final_contigs.taxonomy.tsv",
        gff=results+"/annotation/{assembly}/final_contigs.gff"
    output:
        tax=expand(results+"/annotation/{{assembly}}/taxonomy/orfs.{db}.taxonomy.tsv",
                    db=config["taxonomy"]["database"])
    script:
        "../scripts/taxonomy_utils.py"
