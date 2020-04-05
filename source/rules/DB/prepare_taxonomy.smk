localrules:
    krona_taxonomy,
    tango_download,
    tango_download_taxonomy,
    tango_format_nr,
    tango_format_uniref

rule tango_download_taxonomy:
    output:
        opj(config["resource_path"],"taxonomy","taxonomy.sqlite"),
        opj(config["resource_path"],"taxonomy","taxdump.tar.gz"),
        opj(config["resource_path"],"taxonomy","nodes.dmp"),
        opj(config["resource_path"],"taxonomy","names.dmp"),
        opj(config["resource_path"],"taxonomy","taxonomy.sqlite.traverse.pkl")
    log:
        opj(config["resource_path"],"taxonomy","tango.log")
    params:
        taxdir=opj(config["resource_path"],"taxonomy")
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango download taxonomy -t {params.taxdir} >{log} 2>&1
        """

rule krona_taxonomy:
    output:
        opj(config["resource_path"],"krona","taxonomy.tab")
    log:
        opj(config["resource_path"],"krona","taxonomy.log")
    params:
        taxdir=opj(config["resource_path"],"krona")
    conda:
        "../../../envs/krona.yml"
    shell:
        """
        ktUpdateTaxonomy.sh {params.taxdir} >{log} 2>&1
        """

rule tango_download:
    output:
        fasta = temp(opj(config["resource_path"], "{db}", "{db}.fasta.gz"))
    log:
        opj(config["resource_path"], "{db}", "tango_download.log")
    params:
        dldir = opj(config["resource_path"], "{db}"),
        taxdir = opj(config["resource_path"], "taxonomy")
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango download {wildcards.db} --tmpdir $TMPDIR \
            -d {params.dldir} -t {params.taxdir} --skip_idmap >{log} 2>&1
        """

rule tango_download_nr_idmap:
    output:
        idmap = opj(config["resource_path"], "nr", "prot.accession2taxid.gz")
    log:
        opj(config["resource_path"], "nr", "tango_download_idmap.log")
    params:
        dldir = opj(config["resource_path"], "nr")
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango download idmap -d {params.dldir} > {log} 2>&1
        """

rule tango_format_uniref:
    input:
        fasta = opj(config["resource_path"], "{db}", "{db}.fasta.gz")
    output:
        fasta = opj(config["resource_path"], "{db}", "{db}.reformat.fasta.gz"),
        idmap = opj(config["resource_path"], "{db}", "prot.accession2taxid.gz")
    log:
        opj(config["resource_path"], "{db}", "tango_format.log")
    params:
        tmpdir = config["scratch_path"]
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango format \
            -m {output.idmap} --tmpdir {params.tmpdir} \
            {input.fasta} {output.fasta} > {log} 2>&1
        """

rule tango_format_nr:
    input:
        fasta = opj(config["resource_path"], "nr", "nr.fasta.gz")
    output:
        fasta = opj(config["resource_path"], "nr", "nr.reformat.fasta.gz")
    log:
        opj(config["resource_path"], "nr", "tango_format.log")
    params:
        tmpdir = config["scratch_path"]
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango format \
            --tmpdir {params.tmpdir} 
            {input.fasta} {output.fasta} > {log} 2>&1
        """

rule tango_update:
    input:
        idmap = opj(config["resource_path"], "{db}", "prot.accession2taxid.gz")
    output:
        flag = touch(opj(config["resource_path"], "{db}", "updated"))
    run:
        import os
        if os.path.exists(opj(config["resource_path"], wildcards.db, "idmap.tsv.gz")):
            idfile = opj(config["resource_path"], wildcards.db, "idmap.tsv.gz")
            newfile = "{}.new.gz".format((input.idmap).rstrip(".gz"))
            shell("tango update {input.idmap} {idfile} {newfile}")
            shell("mv {newfile} {input.idmap}")

rule tango_build:
    input:
        flag = opj(config["resource_path"], "{db}", "updated"),
        fasta = opj(config["resource_path"], "{db}", "{db}.reformat.fasta.gz"),
        nodes = opj(config["resource_path"],"taxonomy","nodes.dmp"),
        idmap = opj(config["resource_path"], "{db}", "prot.accession2taxid.gz")
    output:
        opj(config["resource_path"],"{db}","diamond.dmnd")
    log:
        opj(config["resource_path"], "{db}", "diamond.dmnd")
    threads: config["diamond_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../../../envs/tango.yml"
    shell:
         """
         tango build \
            -d {output[0]} \
            -p {threads} {input.fasta} {input.idmap} {input.nodes} 
         """

