import pandas as pd
from scripts.common import annotation_input

localrules:
    annotate,
    download_rfams,
    press_rfams,
    download_pfam,
    download_pfam_info,
    press_pfam,
    parse_pfam,
    download_eggnog,
    get_kegg_info,
    parse_emapper,
    download_rgi_data,
    parse_rgi

##### annotation master rule #####

rule annotate:
    input:
        annotation_input(config, assemblies)

##### gene calling #####

rule prodigal:
    """
    Runs the prodigal gene caller in metagenomic mode on the assembled contigs
    """
    input:
        results+"/assembly/{assembly}/final_contigs.fa"
    output:
        genes=results+"/annotation/{assembly}/final_contigs.ffn",
        faa=results+"/annotation/{assembly}/final_contigs.faa",
        gff=results+"/annotation/{assembly}/final_contigs.gff"
    log:
        results+"/annotation/{assembly}/prodigal.log"
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*2
    conda:
        "../envs/annotation.yml"
    shell:
        """
        prodigal -i {input} -d {output.genes} -a {output.faa} -o {output.gff} \
            -f gff -p meta 2>{log}
        """

rule trnascan:
    input:
        results+"/assembly/{assembly}/final_contigs.fa"
    output:
        file=results+"/annotation/{assembly}/tRNA.out",
        bed=results+"/annotation/{assembly}/tRNA.bed",
        fasta=results+"/annotation/{assembly}/tRNA.fasta"
    log:
        results+"/annotation/{assembly}/tRNA.log"
    threads: 4
    conda:
        "../envs/annotation.yml"
    shell:
        """
        tRNAscan-SE -G -b {output.bed} -o {output.file} -a {output.fasta} \
            --thread {threads} {input} >{log} 2>&1
        """

rule download_rfams:
    output:
        tar=temp("resources/infernal/Rfam.tar.gz"),
        cm="resources/infernal/Rfam.rRNA.cm",
        readme="resources/infernal/README",
        version="resources/infernal/Rfam.version",
        clanin="resources/infernal/Rfam.clanin"
    log:
        "resources/infernal/download.log"
    params:
        url="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT",
        rfams=" ".join(["RF00001.cm", "RF00002.cm", "RF00177.cm",
               "RF01959.cm", "RF01960.cm", "RF02540.cm", "RF02541.cm",
               "RF02542.cm", "RF02543.cm", "RF02545.cm", "RF02546.cm",
               "RF02547.cm", "RF02554.cm", "RF02555.cm"]),
        dir=lambda w, output: os.path.dirname(output.cm)
    shell:
        """
        # Get the Rfam tarball
        curl -o {output.tar} {params.url}/Rfam.tar.gz > {log} 2>&1
        # Extract only rfams of interest
        tar -C {params.dir} -zxf {output.tar} {params.rfams}
        cat {params.dir}/*.cm > {output.cm}

        # Get release
        curl -o {output.readme} {params.url}/README 2>/dev/null
        grep -m 1 Release {output.readme} > {output.version}

        # Get clans
        curl {params.url}/Rfam.clanin 2>/dev/null | egrep -w \
            "CL0011[123]" > {output.clanin}
        """

rule press_rfams:
    input:
        "resources/infernal/Rfam.rRNA.cm"
    output:
        expand("resources/infernal/Rfam.rRNA.cm.i1{suffix}",
               suffix=["m", "i", "f", "p"])
    log:
        "resources/infernal/cmpress.log"
    conda:
        "../envs/annotation.yml"
    shell:
        """
        cmpress {input} > {log} 2>&1
        """

rule infernal:
    input:
        fastafile=results+"/assembly/{assembly}/final_contigs.fa",
        db=expand("resources/infernal/Rfam.rRNA.cm.i1{suffix}",
                  suffix=["m", "i", "f", "p"]),
        cl="resources/infernal/Rfam.clanin"
    output:
        results+"/annotation/{assembly}/final_contigs.cmscan"
    log:
        results+"/annotation/{assembly}/infernal.log"
    params:
        db="resources/infernal/Rfam.rRNA.cm"
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../envs/annotation.yml"
    shell:
        """
        cmscan --cpu {threads} --oskip --rfam --cut_ga --nohmmonly \
            --tblout {output} --fmt 2 --clanin {input.cl} {params.db} \
            {input.fastafile} > /dev/null 2>{log}
        """

##### protein families #####

rule download_pfam:
    output:
        hmmfile="resources/pfam/Pfam-A.hmm",
        datfile="resources/pfam/Pfam-A.hmm.dat",
        versionfile="resources/pfam/Pfam-A.version"
    log:
        "resources/pfam/download.log"
    params:
        ftp="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
    shell:
        """
        curl -s -L -o {output.hmmfile}.gz {params.ftp}/Pfam-A.hmm.gz
        curl -s -L -o {output.datfile}.gz {params.ftp}/Pfam-A.hmm.dat.gz
        curl -s -L -o {output.versionfile}.gz {params.ftp}/Pfam.version.gz

        gunzip {output.hmmfile}.gz
        gunzip {output.datfile}.gz
        gunzip {output.versionfile}.gz
        """

rule download_pfam_info:
    output:
        clanfile="resources/pfam/clan.txt",
        info="resources/pfam/Pfam-A.clans.tsv"
    log:
        "resources/pfam/info.log"
    params:
        ftp="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
    shell:
        """
        curl -s -L -o {output.clanfile}.gz {params.ftp}/database_files/clan.txt.gz
        curl -s -L -o {output.info}.gz {params.ftp}/Pfam-A.clans.tsv.gz

        gunzip {output.clanfile}.gz
        gunzip {output.info}.gz
        """

rule press_pfam:
    input:
        hmmfile="resources/pfam/Pfam-A.hmm"
    output:
        expand("resources/pfam/Pfam-A.hmm.h3{suffix}",
               suffix=["f", "i", "m", "p"])
    log:
        "resources/pfam/hmmpress.log"
    conda:
        "../envs/annotation.yml"
    shell:
        """
        hmmpress {input.hmmfile} > {log} 2>&1
        """

rule pfam_scan:
    input:
        results+"/annotation/{assembly}/final_contigs.faa",
        expand("resources/pfam/Pfam-A.hmm.h3{suffix}",
               suffix=["f", "i", "m", "p"])
    output:
        results+"/annotation/{assembly}/{assembly}.pfam.out"
    log:
        results+"/annotation/{assembly}/{assembly}.pfam.log"
    conda:
        "../envs/annotation.yml"
    params:
        dir="resources/pfam",
        tmp_out=temppath+"/{assembly}.pfam.out"
    threads: 2
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    shell:
        """
        pfam_scan.pl -fasta {input[0]} -dir {params.dir} -cpu {threads} \
            -outfile {params.tmp_out} >{log} 2>&1
        mv {params.tmp_out} {output[0]}
        """

rule parse_pfam:
    input:
        results+"/annotation/{assembly}/{assembly}.pfam.out",
        "resources/pfam/clan.txt",
        "resources/pfam/Pfam-A.clans.tsv"
    output:
        results+"/annotation/{assembly}/pfam.parsed.tsv"
    script:
        "../scripts/annotation_utils.py"

##### eggnog-mapper #####

rule download_eggnog:
    output:
        db = "resources/eggnog-mapper/eggnog.db",
        dmnd = "resources/eggnog-mapper/eggnog_proteins.dmnd",
        version = "resources/eggnog-mapper/eggnog.version"
    log:
        "resources/eggnog-mapper/download.log"
    conda:
        "../envs/annotation.yml"
    params:
        data_dir=lambda wildcards, output: os.path.dirname(output.db)
    shell:
        """
        download_eggnog_data.py --data_dir {params.data_dir} -y > {log} 2>&1
        egrep -o "emapperdb-[0-9].[0-9].[0-9]" {log} > {output.version}
        """

rule get_kegg_info:
    #TODO: Check which files are needed with new eggnog-mapper version
    output:
        kos = "resources/kegg/kegg_kos.tsv",
        mods = "resources/kegg/kegg_modules.tsv",
        pwys = "resources/kegg/kegg_pathways.tsv",
        enzs = touch("resources/kegg/kegg_enzymes.tsv")
    log:
        "resources/kegg/download.log"
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        src="workflow/scripts/eggnog-parser.py"
    shell:
        """
        python {params.src} download {params.outdir} > {log} 2>&1
        """

rule emapper_homology_search:
    input:
        results+"/annotation/{assembly}/final_contigs.faa",
        "resources/eggnog-mapper/eggnog.db",
        "resources/eggnog-mapper/eggnog_proteins.dmnd",
    output:
        results+"/annotation/{assembly}/{assembly}.emapper.seed_orthologs"
    params:
        resource_dir=lambda wildcards, input: os.path.dirname(input[1]),
        out="{assembly}",
        tmpdir=temppath+"/{assembly}-eggnog",
        tmp_out=temppath+"/{assembly}-eggnog/{assembly}",
        flags="-m diamond --no_annot --no_file_comments"
    log:
        results+"/annotation/{assembly}/{assembly}.emapper.seed_orthologs.log"
    conda:
        "../envs/annotation.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        mkdir -p {params.tmpdir}
        emapper.py {params.flags} --cpu {threads} --temp_dir {params.tmpdir} \
        -i {input[0]} -o {params.out} --output_dir {params.tmpdir} \
            --data_dir {params.resource_dir} >{log} 2>&1
        mv {params.tmp_out}.emapper.seed_orthologs {output[0]}
        rm -rf {params.tmpdir}
        """

if config["runOnUppMax"]:
    rule emapper_annotate_hits_uppmax:
        """Copy EGGNOG db into memory before running annotations"""
        input:
            results+"/annotation/{assembly}/{assembly}.emapper.seed_orthologs",
            "resources/eggnog-mapper/eggnog.db",
            "resources/eggnog-mapper/eggnog_proteins.dmnd"
        output:
            results+"/annotation/{assembly}/{assembly}.emapper.annotations"
        params:
            resource_dir=lambda wildcards, input: os.path.dirname(input[1]),
            tmpdir=temppath+"/{assembly}-eggnog",
            out=results+"/annotation/{assembly}/{assembly}",
            flags="--no_file_comments"
        log:
            results+"/annotation/{assembly}/{assembly}.emapper.annotations.log"
        conda:
            "../envs/annotation.yml"
        message: "Annotating hits table for {wildcards.assembly}"
        threads: 10
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60
        shell:
            """
            if [ -z ${{SLURM_JOB_ID+x}} ]; then SLURM_JOB_ID="emapper_annotate_hits_uppmax"; fi
            #Copy eggnog.db
            mkdir -p /dev/shm/$SLURM_JOB_ID
            cp {params.resource_dir}/eggnog.db {params.resource_dir}/eggnog_proteins.dmnd /dev/shm/$SLURM_JOB_ID

            emapper.py {params.flags} --cpu {threads} -o {params.out} \
                --annotate_hits_table {input[0]} --usemem \
                --data_dir /dev/shm/$SLURM_JOB_ID >{log} 2>&1
            rm -rf /dev/shm/$SLURM_JOB_ID
            """
else:
    rule emapper_annotate_hits:
        input:
            results+"/annotation/{assembly}/{assembly}.emapper.seed_orthologs"
        output:
            results+"/annotation/{assembly}/{assembly}.emapper.annotations"
        log:
            results+"/annotation/{assembly}/{assembly}.emapper.annotations.log"
        params:
            resource_dir="resources/eggnog-mapper",
            tmpdir=temppath+"/{assembly}-eggnog",
            out=results+"/annotation/{assembly}/{assembly}",
            flags="--no_file_comments"
        conda:
            "../envs/annotation.yml"
        threads: 10
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60
        shell:
            """
            emapper.py {params.flags} --cpu {threads} -o {params.out} \
                --annotate_hits_table {input[0]} --usemem \
                --data_dir {params.resource_dir} >{log} 2>&1
            """

rule parse_emapper:
    input:
        annotations = results+"/annotation/{assembly}/{assembly}.emapper.annotations",
        info = "resources/kegg/kegg_{db}.tsv"
    wildcard_constraints:
        db="enzymes|pathways|kos|modules"
    output:
        results+"/annotation/{assembly}/{db}.parsed.tsv"
    script:
        "../scripts/annotation_utils.py"

##### resistance gene identifier #####

rule download_rgi_data:
    output:
        json="resources/card/card.json",
        version="resources/card/card.version"
    log:
        "resources/card/log"
    params:
        tar="resources/card/data.tar.gz",
        dir=lambda w, output: os.path.dirname(output.json)
    shell:
         """
         curl -L -o {params.tar} \
            https://card.mcmaster.ca/latest/data >{log} 2>&1
         tar -C {params.dir} -xf {params.tar} ./card.json
         # Store download date in versionfile
         date > {output.version}
         rm {params.tar}
         """

rule rgi:
    input:
        faa=results+"/annotation/{assembly}/final_contigs.faa",
        db="resources/card/card.json"
    output:
        json=results+"/annotation/{assembly}/rgi.out.json",
        txt=results+"/annotation/{assembly}/rgi.out.txt"
    log:
        results+"/annotation/{assembly}/rgi.log"
    params:
        out=results+"/annotation/{assembly}/rgi.out",
        settings="-a diamond --local --clean --input_type protein",
        faa=temppath+"/{assembly}.rgi/final_contig.faa",
        tmpdir=temppath+"/{assembly}.rgi"
    shadow: "minimal"
    conda:
        "../envs/rgi.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    shell:
        """
        mkdir -p {params.tmpdir}
        rgi load -i {input.db} --local > {log} 2>&1
        sed 's/*//g' {input.faa} > {params.faa}
        rgi main -i {params.faa} -o {params.out} \
            -n {threads} {params.settings} >>{log} 2>>{log}
        rm -r {params.tmpdir}
        """

rule parse_rgi:
    input:
        txt=results+"/annotation/{assembly}/rgi.out.txt"
    output:
        tsv=results+"/annotation/{assembly}/rgi.parsed.tsv"
    script:
        "../scripts/annotation_utils.py"
