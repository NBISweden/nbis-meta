from scripts.common import cami_dataset, cami_gold_urls

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

rule download_cami:
    output:
        fq = "data/cami/{c}_S00{s}__insert_{l}_reads_anonymous.fq.gz",
        fa = opj(config["paths"]["results"], "assembly", "{c}_S00{s}__insert_{l}.GOLD", "final_contigs.fa"),
        log = touch("results/assembly/{c}_S00{s}__insert_{l}.GOLD/log")
    log:
        "data/cami/{c}_S00{s}__insert_{l}_reads_anonymous.log"
    params:
        dataset = lambda wildcards, output: cami_dataset(output[0]),
        base_fq = lambda wildcards, output: os.path.basename(output.fq),
        base_fa = lambda wildcards, output: os.path.basename(os.path.dirname(output.fa)).replace(".GOLD", "_gsa_anonymous.fasta.gz"),
        tmp_fq = "$TMPDIR/R{c}_S00{s}__insert_{l}_reads_anonymous.fq.gz",
        tmp_fa = "$TMPDIR/R{c}_S00{s}__insert_{l}_gsa_anonymous.fasta.gz",
        base_url = "https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1"
    shell:
        """
        curl -L -o {params.tmp_fq} {params.base_url}/{params.dataset}/{params.base_fq} > {log} 2>&1
        curl -L -o {params.tmp_fa} {params.base_url}/{params.dataset}/{params.base_fa} > {log} 2>&1
        gunzip -c {params.tmp_fa} > {output.fa}
        mv {params.tmp_fq} {output.fq}
        """

rule download_cami_gold_asms:
    output:
        opj(config["paths"]["results"], "assembly", "{camidataset}.GOLD", "final_contigs.fa"),
        touch(opj(config["paths"]["results"], "assembly", "{camidataset}.GOLD", "log"))
    log:
        opj(config["paths"]["results"], "assembly", "{camidataset}.GOLD", "cami.log")
    params:
        tmpdir = "$TMPDIR/{camidataset}.GOLD",
        tmp = "$TMPDIR/{camidataset}.GOLD/final_contigs.fa.gz",
        url = lambda wildcards: cami_gold_urls(wildcards)
    shell:
        """
        mkdir -p {params.tmpdir}
        curl -L -o {params.tmp} {params.url} > {log} 2>&1
        gunzip -c {params.tmp} > {output[0]}
        rm -r {params.tmpdir}        
        """

ruleorder: download_cami > download_cami_gold_asms > megahit

rule deinterleave_cami_data:
    input:
        "data/cami/{c}_S00{s}__insert_{l}_reads_anonymous.fq.gz"
    output:
        "data/cami/{c}_S00{s}__insert_{l}_reads_anonymous_R{i}.fastq.gz"
    conda:
        "../envs/examples.yml"
    params:
        tmp = "$TMPDIR/{c}_S00{s}__insert_{l}_reads_anonymous_R{i}.fastq.gz"
    shell:
        """
        seqtk seq -{wildcards.i} {input} | \
            sed "s/^@\([0-9A-Za-z]\+\)|\([0-9A-Za-z]\+\)|\([0-9A-Za-z]\+\)\/{wildcards.i}/@\\1|\\2|\\3 {wildcards.i}/g" | \
            gzip -c > {params.tmp}
        mv {params.tmp} {output}
        """
