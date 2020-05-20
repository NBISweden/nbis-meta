localrules:
    tango_download,
    tango_download_taxonomy,
    tango_format_nr,
    tango_format_uniref,
    download_sourmash_db,
    tango_assign_orfs, 
    download_sourmash_db,
    sourmash_compute, 
    merge_tango_sourmash

##### taxonomy master rule #####
rule taxonomy:
    input:
        expand(opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
                   "orfs.{db}.taxonomy.tsv"),
               group=assemblies.keys(), db=config["taxonomy"]["database"])

##### tango #####

rule tango_download_taxonomy:
    output:
        sqlite=opj("resources", "taxonomy", "taxonomy.sqlite"),
        taxdump=opj("resources", "taxonomy", "taxdump.tar.gz"),
        nodes=opj("resources", "taxonomy", "nodes.dmp"),
        names=opj("resources", "taxonomy", "names.dmp"),
        pkl=opj("resources", "taxonomy",
                  "taxonomy.sqlite.traverse.pkl")
    log:
        opj("resources", "taxonomy", "tango.log")
    params:
        taxdir=lambda wildcards, output: os.path.dirname(output.sqlite)
    conda:
        "../envs/tango.yml"
    shell:
        """
        tango download taxonomy -t {params.taxdir} >{log} 2>&1
        """

rule tango_download:
    output:
        fasta=temp(opj("resources", "{db}", "{db}.fasta.gz"))
    log:
        opj("resources", "{db}", "tango_download.log")
    params:
        dldir=lambda wildcards, output: os.path.dirname(output.fasta),
        tmpdir="$TMPDIR"
    conda:
        "../envs/tango.yml"
    shell:
        """
        tango download {wildcards.db} --tmpdir {params.tmpdir} \
            -d {params.dldir} --skip_idmap >{log} 2>{log}
        """

rule tango_download_nr_idmap:
    output:
        idmap=opj("resources", "nr", "prot.accession2taxid.gz")
    log:
        opj("resources", "nr", "tango_download_idmap.log")
    params:
        dldir=lambda wildcards, output: os.path.dirname(output.idmap)
    conda:
        "../envs/tango.yml"
    shell:
        """
        tango download idmap -d {params.dldir} > {log} 2>&1
        """

rule tango_format_uniref:
    input:
        fasta=opj("resources", "{db}", "{db}.fasta.gz")
    output:
        fasta=opj("resources", "{db}", "{db}.reformat.fasta.gz"),
        idmap=opj("resources", "{db}", "prot.accession2taxid.gz")
    log:
        opj("resources", "{db}", "tango_format.log")
    params:
        tmpdir=config["paths"]["temp"]
    conda:
        "../envs/tango.yml"
    shell:
        """
        tango format -m {output.idmap} --tmpdir {params.tmpdir} \
            {input.fasta} {output.fasta} > {log} 2>&1
        """

rule tango_format_nr:
    input:
        fasta=opj("resources", "nr", "nr.fasta.gz")
    output:
        fasta=opj("resources", "nr", "nr.reformat.fasta.gz")
    log:
        opj("resources", "nr", "tango_format.log")
    params:
        tmpdir=config["paths"]["temp"]
    conda:
        "../envs/tango.yml"
    shell:
        """
        tango format --tmpdir {params.tmpdir} {input.fasta} \
            {output.fasta} > {log} 2>&1
        """

rule tango_update:
    input:
        idmap=opj("resources", "{db}", "prot.accession2taxid.gz")
    output:
        idmap=opj("resources", "{db}", "prot.accession2taxid.update.gz")
    log:
        opj("resources", "{db}", "tango_update.log")
    conda:
        "../envs/tango.yml"
    params:
        dir=lambda wildcards, output: os.path.dirname(output.idmap)
    shell:
        """
        # If an idmap file is available, use it to create an updated idmap file
        if [ -e {params.dir}/idmap.tsv.gz ] ; then
            tango update {input.idmap} {params.dir}/idmap.tsg.gz \
                {output.idmap} > {log} 2>&1
        # Otherwise, just create a symlink
        else
            cd {params.dir}
            ln -s $(basename {input.idmap}) $(basename {output.idmap})
        fi            
        """

rule tango_build:
    input:
        fasta=opj("resources", "{db}", "{db}.reformat.fasta.gz"),
        nodes=opj("resources", "taxonomy", "nodes.dmp"),
        idmap=opj("resources", "{db}", "prot.accession2taxid.update.gz")
    output:
        opj("resources", "{db}", "diamond.dmnd")
    log:
        opj("resources", "{db}", "diamond.log")
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/tango.yml"
    shell:
         """
         tango build -d {output} -p {threads} {input.fasta} \
            {input.idmap} {input.nodes} >{log} 2>&1 
         """

rule tango_search:
    input:
        db=opj("resources", config["taxonomy"]["database"], "diamond.dmnd"),
        fasta=opj(config["paths"]["results"], "assembly", "{group}",
                  "final_contigs.fa")
    output:
        opj(config["paths"]["results"], "annotation", "{group}",
            "final_contigs.{db}.tsv.gz".format(db=config["taxonomy"]["database"]))
    log:
        opj(config["paths"]["results"], "annotation", "{group}", "tango_search.log")
    params:
        tmpdir=config["paths"]["temp"],
        min_len=config["taxonomy"]["min_len"],
        settings=config["taxonomy"]["search_params"]
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/tango.yml"
    shell:
        """
        tango search {params.settings} -p {threads} \
            --tmpdir {params.tmpdir} -l {params.min_len} \
            {input.fasta} {input.db} {output} >{log} 2>&1
        """

rule tango_assign:
    input:
        tsv=opj(config["paths"]["results"], "annotation", "{group}",
            "final_contigs.{db}.tsv.gz".format(db=config["taxonomy"]["database"])),
        sql=ancient(opj("resources", "taxonomy", "taxonomy.sqlite"))
    output:
        opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
            "tango.{db}.taxonomy.tsv".format(db=config["taxonomy"]["database"]))
    log:
        opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
            "tango_assign.log")
    params:
        taxonomy_ranks=" ".join(config["taxonomy"]["ranks"]),
        taxdir=opj("resources", "taxonomy"),
        settings=config["taxonomy"]["assign_params"]
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../envs/tango.yml"
    shell:
         """
         tango assign {params.settings} -p {threads} -m rank_lca \
            --reportranks {params.taxonomy_ranks} -t {params.taxdir} \
            {input.tsv} {output} > {log} 2>&1
         """

##### sourmash #####

rule download_sourmash_db:
    output:
        opj("resources", "sourmash", "genbank-k31.lca.json")
    log:
        opj("resources", "sourmash", "download.log")
    params:
        url="https://osf.io/4f8n3/download"
    shell:
        """
        curl -L -o {output}.gz {params.url} > {log} 2>&1
        gunzip {output}.gz
        """

rule sourmash_compute:
    input:
        opj(config["paths"]["results"], "assembly", "{group}", "final_contigs.fa")
    output:
        opj(config["paths"]["results"], "assembly", "{group}", "final_contigs.fa.sig")
    log:
        opj(config["paths"]["results"], "assembly", "{group}", "sourmash_compute.log")
    conda:
        "../envs/sourmash.yml"
    params:
        frac=config["taxonomy"]["sourmash_fraction"]
    shell:
        """
        sourmash compute --singleton --scaled {params.frac} \
            -k 31 -o {output} {input} > {log} 2>&1
        """

rule sourmash_classify:
    input:
        sig=opj(config["paths"]["results"], "assembly", "{group}",
                 "final_contigs.fa.sig"),
        db=opj("resources", "sourmash", "genbank-k31.lca.json")
    output:
        csv=opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
                  "sourmash.taxonomy.csv")
    log:
        opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
            "sourmash.log")
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

rule merge_tango_sourmash:
    input:
        smash=opj(config["paths"]["results"], "annotation", "{group}",
                    "taxonomy", "sourmash.taxonomy.csv"),
        tango=opj(config["paths"]["results"], "annotation", "{group}",
                    "taxonomy", "tango.{db}.taxonomy.tsv".format(db=config["taxonomy"]["database"]))
    output:
        opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
        "final_contigs.taxonomy.tsv")
    log:
        opj(config["paths"]["results"], "annotation", "{group}", "taxonomy", "merge.log")
    script:
        "../scripts/taxonomy_utils.py"

rule tango_assign_orfs:
    input:
        tax=opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
            "final_contigs.taxonomy.tsv"),
        gff=opj(config["paths"]["results"], "annotation", "{group}",
                "final_contigs.gff")
    output:
        tax=opj(config["paths"]["results"], "annotation", "{group}", "taxonomy",
            "orfs.{db}.taxonomy.tsv".format(db=config["taxonomy"]["database"]))
    script:
        "../scripts/taxonomy_utils.py"
