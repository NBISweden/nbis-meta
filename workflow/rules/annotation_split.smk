scattergather:
    split = config["annotation"]["splits"]

localrules:
    split_fasta,
    pfam_scan_gather

rule split_fasta:
    input:
        "results/annotation/{assembly}/final_contigs.faa"
    output:
        expand("results/annotation/{{assembly}}/splits/split_{n}-of-{N}.faa",
        n = list(range(1, config["annotation"]["splits"]+1)),
        N =  config["annotation"]["splits"])
    params:
        n_files = lambda wildcards: config["annotation"]["splits"],
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    log: "results/annotation/{assembly}/splits/splits.log"
    shell:
        "python workflow/scripts/split_fasta_file.py {input} -n {params.n_files} -o {params.outdir} > {log} 2>&1"

rule pfam_scan_split:
    input:
        "results/annotation/{assembly}/splits/split_{scatteritem}.faa",
        expand("resources/pfam/Pfam-A.hmm.h3{suffix}",
               suffix=["f", "i", "m", "p"])
    output:
        "results/annotation/{assembly}/splits/split_{scatteritem}.pfam.out"
    conda:
        "../envs/annotation.yml"
    params:
        dir="resources/pfam",
        tmp_out=temppath+"/{assembly}.pfam.out"
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    shell:
        """
        pfam_scan.pl -fasta {input[0]} -dir {params.dir} -cpu {threads} \
            -outfile {params.tmp_out} >{log} 2>&1
        mv {params.tmp_out} {output[0]}
        """

rule pfam_scan_gather:
    input:
        gather.split("results/annotation/{{assembly}}/splits/split_{scatteritem}.pfam.out")
    output:
        "results/annotation/{assembly}/{assembly}.pfam.out",
        touch("results/annotation/{assembly}/{assembly}.pfam.gathered")
    shell:
        """
        egrep "^#" {input[0]} > {output[0]}
        cat {input} | egrep -v "^#" >> {output[0]}
        """
