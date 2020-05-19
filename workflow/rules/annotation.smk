from scripts.common import annotation_input

localrules: 
    parse_ko_annotations, 
    parse_pfam,
    download_rfams,
    press_rfams,
    download_pfam,
    press_pfam,
    download_pfam_info,
    download_eggnog,
    get_kegg_info,
    download_rgi_data

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
        opj(config["results_path"], "assembly", "{group}", "final_contigs.fa")
    output:
        genes=opj(config["results_path"], "annotation", "{group}", "final_contigs.ffn"),
        faa=opj(config["results_path"], "annotation", "{group}", "final_contigs.faa"),
        gff=opj(config["results_path"], "annotation", "{group}", "final_contigs.gff")
    log:
        opj(config["results_path"], "annotation", "{group}", "prodigal.log")
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
        opj(config["results_path"], "assembly", "{group}", "final_contigs.fa")
    output:
        file=opj(config["results_path"], "annotation", "{group}", "tRNA.out"),
        bed=opj(config["results_path"], "annotation", "{group}", "tRNA.bed"),
        fasta=opj(config["results_path"], "annotation", "{group}", "tRNA.fasta")
    log:
        opj(config["results_path"], "annotation", "{group}", "tRNA.log")
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
        tar=temp(opj(config["resource_path"], "infernal", "Rfam.tar.gz")),
        cm=opj(config["resource_path"], "infernal", "Rfam.rRNA.cm"),
        readme=opj(config["resource_path"], "infernal", "README"),
        version=opj(config["resource_path"], "infernal", "Rfam.version"),
        clanin=opj(config["resource_path"], "infernal", "Rfam.clanin")
    log:
        opj(config["resource_path"], "infernal", "download.log")
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
        opj(config["resource_path"], "infernal", "Rfam.rRNA.cm")
    output:
        expand(opj(config["resource_path"], "infernal", "Rfam.rRNA.cm.i1{suffix}"),
               suffix=["m", "i", "f", "p"])
    log:
        opj(config["resource_path"], "infernal", "cmpress.log")
    conda:
        "../envs/annotation.yml"
    shell:
        """
        cmpress {input} > {log} 2>&1
        """

rule infernal:
    input:
        fastafile=opj(config["results_path"], "assembly", "{group}",
                      "final_contigs.fa"),
        db=expand(opj(config["resource_path"], "infernal",
                      "Rfam.rRNA.cm.i1{suffix}"),
                  suffix=["m", "i", "f", "p"]),
        cl=opj(config["resource_path"], "infernal", "Rfam.clanin")
    output:
        opj(config["results_path"], "annotation", "{group}",
            "final_contigs.cmscan")
    log:
        opj(config["results_path"], "annotation", "{group}",
            "infernal.log")
    params:
        db=opj(config["resource_path"], "infernal", "Rfam.rRNA.cm")
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
        hmmfile=opj(config["resource_path"], "pfam", "Pfam-A.hmm"),
        datfile=opj(config["resource_path"], "pfam", "Pfam-A.hmm.dat"),
        versionfile=opj(config["resource_path"], "pfam", "Pfam-A.version"),
    log:
        opj(config["resource_path"], "pfam", "download.log")
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
        clanfile=opj(config["resource_path"], "pfam", "clan.txt"),
        info=opj(config["resource_path"], "pfam", "Pfam-A.clans.tsv")
    log:
        opj(config["resource_path"], "pfam", "info.log")
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
        hmmfile=opj(config["resource_path"], "pfam", "Pfam-A.hmm")
    output:
        expand(opj(config["resource_path"], "pfam", "Pfam-A.hmm.h3{suffix}"),
               suffix=["f", "i", "m", "p"])
    log:
        opj(config["resource_path"], "pfam", "hmmpress.log")
    conda:
        "../envs/annotation.yml"
    shell:
        """
        hmmpress {input.hmmfile} > {log} 2>&1
        """

rule pfam_scan:
    input:
        opj(config["results_path"], "annotation", "{group}", "final_contigs.faa"),
        expand(opj(config["resource_path"], "pfam", "Pfam-A.hmm.h3{suffix}"),
               suffix=["f", "i", "m", "p"])
    output:
        opj(config["results_path"], "annotation", "{group}", "{group}.pfam.out")
    log:
        opj(config["results_path"], "annotation", "{group}", "{group}.pfam.log")
    conda:
        "../envs/annotation.yml"
    params:
        dir=opj(config["resource_path"], "pfam"),
        tmp_out=opj(os.path.expandvars(config["scratch_path"]),
                    "{group}.pfam.out")
    threads: 2
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        pfam_scan.pl -fasta {input[0]} -dir {params.dir} -cpu {threads} \
            -outfile {params.tmp_out} >{log} 2>&1
        mv {params.tmp_out} {output[0]}
        """

rule parse_pfam:
    input:
        opj(config["results_path"], "annotation", "{group}", "{group}.pfam.out"),
        opj(config["resource_path"], "pfam", "clan.txt"),
        opj(config["resource_path"], "pfam", "Pfam-A.clans.tsv")
    output:
        opj(config["results_path"], "annotation", "{group}", "pfam.parsed.tsv")
    script:
        "../scripts/annotation_utils.py"

##### eggnog-mapper #####

rule download_eggnog:
    output:
        db=opj(config["resource_path"],"eggnog-mapper","eggnog.db"),
        version=opj(config["resource_path"],"eggnog-mapper","eggnog.version")
    log:
        opj(config["resource_path"],"eggnog-mapper","download.log")
    conda:
        "../envs/annotation.yml"
    params:
        dbs="none",
        data_dir=lambda wildcards, output: os.path.dirname(output.db)
    shell:
        """
        download_eggnog_data.py --data_dir {params.data_dir} -y > {log} 2>&1
        egrep -o "emapperdb-[0-9].[0-9].[0-9]" {log} > {output.version}
        """

rule get_kegg_info:
    #TODO: Check which files are needed with new eggnog-mapper version
    output:
        expand(opj(config["resource_path"], "kegg", "{f}"),
               f=["kegg_ec2pathways.tsv", "kegg_ko2ec.tsv",
                  "kegg_ko2pathways.tsv", "kegg_kos.tsv", "kegg_modules.tsv",
                  "kegg_pathways.tsv"])
    log:
        opj(config["resource_path"], "kegg", "download.log")
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        src="../scripts/eggnog-parser.py"
    shell:
        """
        python {params.src} download {params.outdir} > {log} 2>&1
        """

rule emapper_homology_search:
    input:
        opj(config["results_path"], "annotation", "{group}",
            "final_contigs.faa"),
        opj(config["resource_path"], "eggnog-mapper", "eggnog.db")
    output:
        opj(config["results_path"], "annotation", "{group}",
            "{group}.emapper.seed_orthologs")
    params:
        resource_dir=opj(config["resource_path"], "eggnog-mapper"),
        out="{group}",
        tmpdir=opj(os.path.expandvars(config["scratch_path"]),
                   "{group}-eggnog"),
        tmp_out=opj(os.path.expandvars(config["scratch_path"]),
                    "{group}-eggnog", "{group}"),
        flags="-m diamond --no_annot --no_file_comments"
    log:
        opj(config["results_path"], "annotation", "{group}",
            "{group}.emapper.seed_orthologs.log")
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

if config["runOnUppMax"] == "yes":
    rule emapper_annotate_hits_uppmax:
        """Copy EGGNOG db into memory before running annotations"""
        input:
            opj(config["results_path"], "annotation", "{group}",
                "{group}.emapper.seed_orthologs")
        output:
            opj(config["results_path"], "annotation", "{group}",
                "{group}.emapper.annotations")
        params:
            resource_dir=opj(config["resource_path"], "eggnog-mapper"),
            tmpdir=opj(os.path.expandvars(config["scratch_path"]),
                       "{group}-eggnog"),
            out=opj(config["results_path"], "annotation", "{group}", "{group}"),
            flags="--no_file_comments"
        log:
            opj(config["results_path"], "annotation", "{group}",
                 "{group}.emapper.annotations.log")
        conda:
            "../envs/annotation.yml"
        message: "Annotating hits table for {wildcards.group}"
        threads: 10
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60
        shell:
            """
            #Copy eggnog.db
            mkdir -p /dev/shm/$SLURM_JOB_ID
            cp {params.resource_dir}/eggnog.db /dev/shm/$SLURM_JOB_ID
            emapper.py {params.flags} --cpu {threads} -o {params.out} \
                --annotate_hits_table {input[0]} --usemem \
                --data_dir /dev/shm/$SLURM_JOB_ID >{log} 2>&1
            rm -rf /dev/shm/$SLURM_JOB_ID
            """
else:
    rule emapper_annotate_hits:
        input:
            opj(config["results_path"], "annotation", "{group}",
                "{group}.emapper.seed_orthologs")
        output:
            opj(config["results_path"], "annotation", "{group}",
                "{group}.emapper.annotations")
        log:
            opj(config["results_path"], "annotation", "{group}",
                 "{group}.emapper.annotations.log")
        params:
            resource_dir=opj(config["resource_path"], "eggnog-mapper"),
            tmpdir=opj(os.path.expandvars(config["scratch_path"]),
                       "{group}-eggnog"),
            out=opj(config["results_path"], "annotation", "{group}", "{group}"),
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

rule parse_ko_annotations:
    input:
        annotations=opj(config["results_path"], "annotation", "{group}", "{group}.emapper.annotations"),
        ko2ec=opj(config["resource_path"], "kegg", "kegg_ko2ec.tsv"),
        ko2path=opj(config["resource_path"], "kegg", "kegg_ko2pathways.tsv"),
        #ko2module=opj(config["resource_path"], "kegg", "kegg_ko2modules.tsv"),
        kos=opj(config["resource_path"], "kegg", "kegg_kos.tsv"),
        modules=opj(config["resource_path"], "kegg", "kegg_modules.tsv"),
        pathways=opj(config["resource_path"], "kegg", "kegg_pathways.tsv")
    output:
        expand(opj(config["results_path"], "annotation", "{{group}}", "{db}.parsed.tsv"),
            db=["enzymes", "pathways", "modules", "kos"])
    log:
        opj(config["results_path"], "annotation", "{{group}}", "eggnog-parser.log")
    params:
        outbase=opj(config["results_path"], "annotation", "{group}"),
        resource_dir=opj(config["resource_path"], "kegg"),
        src="../scripts/eggnog-parser.py"
    shell:
        """
        python {params.src} parse {params.resource_dir} {input.annotations} \
            {params.outbase} > {log} 2>&1
        """

##### resistance gene identifier #####

rule download_rgi_data:
    output:
        json=opj(config["resource_path"], "card", "card.json"),
        version=opj(config["resource_path"], "card", "card.version")
    log:
        opj(config["resource_path"], "card", "log")
    params:
        tar=opj(config["resource_path"], "card", "data.tar.gz"),
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
        faa=opj(config["results_path"], "annotation", "{group}", "final_contigs.faa"),
        db=opj(config["resource_path"], "card", "card.json")
    output:
        json=opj(config["results_path"], "annotation", "{group}", "rgi.out.json"),
        txt=opj(config["results_path"], "annotation", "{group}", "rgi.out.txt")
    log:
        opj(config["results_path"], "annotation", "{group}", "rgi.log")
    params:
        out=opj(config["results_path"], "annotation", "{group}", "rgi.out"),
        settings="-a diamond --local --clean --input_type protein"
    conda:
        "../envs/rgi.yml"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60
    shell:
        """
        rgi load --card_json {input.db} --local > {log} 2>&1
        rgi main -i {input.faa} -o {params.out} \
            -n {threads} {params.settings} >>{log} 2>>{log}
        """
