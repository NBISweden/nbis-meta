localrules:
    download_rRNA_database

rule download_rRNA_database:
    output:
        opj(config["resource_path"],"rRNA_databases","{file}.fasta")
    log:
        opj(config["resource_path"],"rRNA_databases","{file}.dl.log")
    params:
        url="https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/{file}.fasta"
    shell:
         """
         curl -L -o {output} {params.url} > {log} 2>&1
         """

rule index_db:
    input:
        fasta=opj(config["resource_path"],"rRNA_databases","{file}")
    output:
        expand(opj(config["resource_path"],"rRNA_databases",
            "{{file}}.{suffix}"), suffix=["bursttrie_0.dat","kmer_0.dat",
                                              "pos_0.dat","stats"])
    log:
        opj(config["resource_path"], "rRNA_databases", "{file}.index.log")
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*5
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        indexdb_rna --ref {input.fasta},{input.fasta} > {log} 2>&1
        """