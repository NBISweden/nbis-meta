##### download example files #####

rule download_synthetic:
    """
    Download pre-made synthetic metagenome from Zenodo
    """
    output:
        R1 = temp("examples/data/synthetic_1.fastq.gz"),
        R2 = temp("examples/data/synthetic_2.fastq.gz")
    log:
        "examples/data/synthetic.log"
    params:
        tar = "examples/data/synthetic.tar.gz",
        url = "https://zenodo.org/record/3737112/files/synthetic.tar.gz?download=1",
        outdir = lambda wildcards, output: os.path.dirname(output.R1)
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
        "examples/data/synthetic_{i}.fastq.gz"
    output:
        "examples/data/{sample}_{s}_R{i}.fastq.gz"
    log:
        "examples/data/{sample}_{s}_R{i}.log"
    conda:
        "../envs/examples.yml"
    params:
        example_dataset_size=config["example_dataset_size"]
    shell:
         """
         seqtk sample -s {wildcards.s} {input} \
            {params.example_dataset_size} | gzip -c > {output}
         """

rule download_cami_data:
    output:
        "data/cami/{dataset}.fq.gz"
    shadow: "minimal"
    shell:
        """
        # Download cami client
        curl -O https://data.cami-challenge.org/camiClient.jar
        java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/{wildcards.dataset} $TMPDIR/ -p fq.gz
        mv $TMPDIR/*.fq.gz {output[0]}
        """

rule deinterleave_cami_data:
    input:
        "data/cami/{dataset}.fq.gz"
    output:
        "data/cami/{dataset}_R{i}.fastq.gz"
    conda:
        "../envs/examples.yml"
    shell:
        """
        seqtk seq -{wildcards.i} {input} | \
            sed 's/^@\([0-9A-Za-z]\+\)|\([0-9A-Za-z]\+\)|\([0-9A-Za-z]\+\)\/{wildcards.i}/@\1|\2|\3 {wildcards.i}/g' | \
            gzip -c > $TMPDIR/{wildcards.dataset}_R{wildcards.i}.fastq.gz
        mv $TMPDIR/{wildcards.dataset}_R{wildcards.i}.fastq.gz {output}
        """
