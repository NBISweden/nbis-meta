##### download example files #####
localrules:
    download_synthetic,
    generate_examples

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
