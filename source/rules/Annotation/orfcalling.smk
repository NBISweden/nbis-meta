localrules:
    infernal_update_gff

rule prodigal:
    """
    Runs the prodigal gene caller in metagenomic mode on the assembled contigs
    """
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        genes=opj(config["results_path"],"annotation","{group}","final_contigs.ffn"),
        faa=opj(config["results_path"],"annotation","{group}","final_contigs.faa"),
        gff=opj(config["results_path"],"annotation","{group}","final_contigs.gff")
    log:
        opj(config["results_path"],"annotation","{group}","prodigal.log")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*2
    conda:
        "../../../envs/annotation.yml"
    shell:
        """
        prodigal \
            -i {input} \
            -d {output.genes} \
            -a {output.faa} \
            -o {output.gff} \
            -f gff \
            -p meta 2>{log}
        """

rule trnascan:
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        file=opj(config["results_path"],"annotation","{group}",
                 "tRNA.out"),
        bed=opj(config["results_path"],"annotation","{group}",
                 "tRNA.bed"),
        fasta=opj(config["results_path"],"annotation","{group}",
                 "tRNA.fasta")
    log:
        opj(config["results_path"],"annotation","{group}",
                 "tRNA.log")
    threads: 4
    conda:
        "../../../envs/annotation.yml"
    shell:
        """
        tRNAscan-SE \
            -G -b {output.bed} -o {output.file} -a {output.fasta} \
            --thread {threads} {input} >{log} 2>&1            
        """


rule infernal:
    input:
        fastafile=opj(config["results_path"],"assembly","{group}",
                      "final_contigs.fa"),
        db=expand(opj(config["infernal_dbpath"],"Rfam.rRNA.cm.i1{suffix}"),
               suffix=["m","i","f","p"]),
        cl=opj(config["infernal_dbpath"],"Rfam.clanin")
    output:
        opj(config["results_path"],"annotation","{group}",
            "final_contigs.cmscan")
    log:
        opj(config["results_path"],"annotation","{group}",
            "infernal.log")
    params:
        db=opj(config["infernal_dbpath"],"Rfam.rRNA.cm")
    threads: config["infernal_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../../../envs/annotation.yml"
    shell:
        """
        cmscan \
            --cpu {threads} \
            --oskip \
            --rfam \
            --cut_ga \
            --nohmmonly \
            --tblout {output} \
            --fmt 2 \
            --clanin {input.cl} \
            {params.db} \
            {input.fastafile} > /dev/null 2>{log}
        """


rule infernal_update_gff:
    input:
        gff=opj(config["results_path"],"annotation","{group}","final_contigs.gff"),
        cmout=opj(config["results_path"],"annotation","{group}","final_contigs.rfam"),
    output:
        gff=opj(config["results_path"],"annotation","{group}","final_contigs.infernal.gff")
    run:
        import gffutils
        lines = parse_cmout(input.cmout)
        features = []
        for line in lines:
            feature = gffutils.feature.feature_from_line(line, strict = False, dialect = gffutils.constants.dialect)
            features.append(feature)
        gffdb = gffutils.create_db(input.gff, "{}.db".format(input.gff),
                                   force=True)
        gffdb.update(features)
        with open(output.gff, 'w') as fh:
            for feature in gffdb.all_features():
                fh.write("{}\n".format(str(feature)))