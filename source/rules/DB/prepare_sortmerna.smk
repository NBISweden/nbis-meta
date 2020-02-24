localrules:
    download_rRNA_database

rule download_rRNA_database:
    output:
        opj(config["resource_path"],"rRNA_databases","{file}.fasta")
    params:
        url="https://raw.githubusercontent.com/biocore/sortmerna/master/data"
            "/rRNA_databases/{file}.fasta"
    shell:
         """
         curl -L -o {output[0]} {params.url}
         """

rule index_db:
    input:
        fasta=opj(config["resource_path"],"rRNA_databases","{file}")
    output:
        expand(opj(config["resource_path"],"rRNA_databases",
            "{{file}}.{suffix}"), suffix=["bursttrie_0.dat","kmer_0.dat",
                                              "pos_0.dat","stats"])
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*5
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        indexdb_rna --ref {input.fasta},{input.fasta}
        """