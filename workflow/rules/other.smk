##### download example files #####
localrules:
    download_synthetic,
    generate_examples,
    download_cami,
    deinterleave_cami_data


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

def cami_dataset(f):
    cami_datasets = {
        'RL_S001__insert_270.fq.gz': 'CAMI_I_LOW',
        'RM1_S001__insert_5000.fq.gz': 'CAMI_I_MEDIUM',
        'RM1_S002__insert_5000.fq.gz': 'CAMI_I_MEDIUM',
        'RM2_S001__insert_270.fq.gz': 'CAMI_I_MEDIUM',
        'RM2_S002__insert_270.fq.gz': 'CAMI_I_MEDIUM',
        'RH1_S001__insert_270.fq.gz': 'CAMI_I_HIGH',
        'RH1_S002__insert_270.fq.gz': 'CAMI_I_HIGH',
        'RH1_S003__insert_270.fq.gz': 'CAMI_I_HIGH',
        'RH1_S004__insert_270.fq.gz': 'CAMI_I_HIGH',
        'RH1_S005__insert_270.fq.gz': 'CAMI_I_HIGH'
    }
    return cami_datasets[os.path.basename(f)]

rule download_cami:
    output:
        "data/cami/R{c}_S00{s}__insert_{l}.fq.gz"
    log:
        "data/cami/R{c}_S00{s}__insert_{l}.log"
    params:
        dataset = lambda wildcards, output: cami_dataset(output[0]),
        base = lambda wildcards, output: os.path.basename(output[0]),
        tmp = "$TMPDIR/R{c}_S00{s}__insert_{l}.fq.gz"
    shell:
        """
        curl -L -o {params.tmp} https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/{params.dataset}/{params.base} > {log} 2>&1
        mv {params.tmp} {output}
        """

rule deinterleave_cami_data:
    input:
        "data/cami/R{c}_S00{s}__insert_{l}.fq.gz"
    output:
        "data/cami/R{c}_S00{s}__insert_{l}_R{i}.fastq.gz"
    conda:
        "../envs/examples.yml"
    params:
        tmp = "$TMPDIR/R{c}1_S00{s}__insert_{l}_R{i}.fastq.gz"
    shell:
        """
        seqtk seq -{wildcards.i} {input} | \
            sed "s/^@\([0-9A-Za-z]\+\)|\([0-9A-Za-z]\+\)|\([0-9A-Za-z]\+\)\/{wildcards.i}/@\\1|\\2|\\3 {wildcards.i}/g" | \
            gzip -c > {params.tmp}
        mv {params.tmp} {output}
        """
