
rule download_synthetic:
    """
    Download pre-made synthetic metagenome from Zenodo
    """
    output:
        R1 = temp("examples/data/synthetic_1.fastq.gz"),
        R2 = temp("examples/data/synthetic_2.fastq.gz")
    params:
        tar = "examples/data/synthetic.tar.gz",
        url = "https://zenodo.org/record/3737112/files/synthetic.tar.gz?download=1"
    conda:
        "../../envs/examples.yml"
    shell:
         """
         curl -L -s -o {params.tar} {params.url}
         tar -C examples/data/ -xf {params.tar}
         rm {params.tar}
         """

rule generate_examples:
    """
    Use seqtk to subsample the synthetic metagenome into examples
    """
    input:
        "examples/data/synthetic_{i}.fastq.gz"
    output:
        "examples/data/{sample}_{s}_R{i}.fastq.gz"
    conda:
        "../../envs/examples.yml"
    shell:
         """
         seqtk sample -s {wildcards.s} {input[0]} 100000 | gzip -c > {output[0]}
         """