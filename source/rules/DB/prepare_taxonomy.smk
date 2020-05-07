localrules:
    krona_taxonomy,
    tango_download,
    tango_download_taxonomy,
    tango_format_nr,
    tango_format_uniref

rule tango_download_taxonomy:
    output:
        sqlite = opj(config["resource_path"], "taxonomy", "taxonomy.sqlite"),
        taxdump = opj(config["resource_path"], "taxonomy", "taxdump.tar.gz"),
        nodes = opj(config["resource_path"], "taxonomy", "nodes.dmp"),
        names = opj(config["resource_path"], "taxonomy", "names.dmp"),
        pkl = opj(config["resource_path"], "taxonomy",
                  "taxonomy.sqlite.traverse.pkl")
    log:
        opj(config["resource_path"],"taxonomy","tango.log")
    params:
        taxdir = lambda w, output: os.path.dirname(output.sqlite)
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango download taxonomy -t {params.taxdir} >{log} 2>&1
        """

rule krona_taxonomy:
    output:
        tab = opj(config["resource_path"],"krona","taxonomy.tab")
    log:
        opj(config["resource_path"],"krona","taxonomy.log")
    params:
        taxdir = lambda w, output: os.path.dirname(output.tab)
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
        dldir = lambda w, output: os.path.dirname(output.fasta),
        tmpdir = "$TMPDIR"
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango download {wildcards.db} --tmpdir {params.tmpdir} \
            -d {params.dldir} --skip_idmap >{log} 2>{log}
        """

rule tango_download_nr_idmap:
    output:
        idmap = opj(config["resource_path"], "nr", "prot.accession2taxid.gz")
    log:
        opj(config["resource_path"], "nr", "tango_download_idmap.log")
    params:
        dldir = lambda w, output: os.path.dirname(output.idmap)
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
            --tmpdir {params.tmpdir} \
            {input.fasta} {output.fasta} > {log} 2>&1
        """

rule tango_update:
    input:
        idmap = opj(config["resource_path"], "{db}", "prot.accession2taxid.gz")
    output:
        idmap = opj(config["resource_path"], "{db}", "prot.accession2taxid.update.gz")
    log:
        opj(config["resource_path"], "{db}", "tango_update.log")
    conda:
        "../../../envs/tango.yml"
    params:
        dir = lambda w, output: os.path.dirname(output.idmap)
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
        fasta = opj(config["resource_path"], "{db}", "{db}.reformat.fasta.gz"),
        nodes = opj(config["resource_path"],"taxonomy","nodes.dmp"),
        idmap = opj(config["resource_path"], "{db}", "prot.accession2taxid.update.gz")
    output:
        opj(config["resource_path"],"{db}","diamond.dmnd")
    log:
        opj(config["resource_path"], "{db}", "diamond.log")
    threads: config["diamond_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../../../envs/tango.yml"
    shell:
         """
         tango build -d {output} -p {threads} {input.fasta} \
            {input.idmap} {input.nodes} >{log} 2>&1 
         """