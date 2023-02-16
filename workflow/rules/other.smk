##### download example files #####
localrules:
    download_synthetic,
    generate_examples,


rule download_synthetic:
    """
    Download pre-made synthetic metagenome from Zenodo
    """
    output:
        R1=temp("examples/data/synthetic_1.fastq.gz"),
        R2=temp("examples/data/synthetic_2.fastq.gz"),
    log:
        "examples/data/synthetic.log",
    params:
        tar="examples/data/synthetic.tar.gz",
        url="https://zenodo.org/record/3737112/files/synthetic.tar.gz?download=1",
        outdir=lambda wildcards, output: os.path.dirname(output.R1),
    shell:
        """
         curl -L -s -o {params.tar} {params.url}
         tar -C {params.outdir} -xf {params.tar}
         rm {params.tar}
         """


rule generate_examples:
    """
    Use seqtk to subsample the synthetic metagenome into examples
    """
    input:
        "examples/data/synthetic_{i}.fastq.gz",
    output:
        "examples/data/{sample}_{s}_R{i}.fastq.gz",
    log:
        "examples/data/{sample}_{s}_R{i}.log",
    conda:
        "../envs/examples.yml"
    envmodules:
        "bioinfo-tools",
        "seqtk",
    params:
        example_dataset_size=config["example_dataset_size"],
    shell:
        """
         seqtk sample -s {wildcards.s} {input} \
            {params.example_dataset_size} | gzip -c > {output}
         """

rule prep_taxonomy:
    output:
        expand("resources/mmseqs/UniRef100{suff}",
            suff=["", "_taxonomy", "_mapping", "_h.dbtype", "_h.index", "_h",
                ".lookup", ".dbtype", ".index", ])
    input:
        fasta=".test/data/uniref100.fasta",
        taxmap=".test/data/taxidmap"
    params:
        tmpdir = "$TMPDIR"
    conda:
        "../envs/mmseqs.yml"
    shell:
        """
        mmseqs createdb {input.fasta} {output[0]}
        mmseqs createtaxdb {output[0]} {params.tmpdir} --tax-mapping-file {input.taxmap}
        """
