from scripts.common import get_all_assembly_files, get_bamfiles

localrules:
    assemble,
    fasta2bed,
    plot_assembly_stats,
    assembly_stats,
    samtools_flagstat

##### master assembly rule #####

rule assemble:
    input:
        expand("{results_path}/report/assembly/{f}.pdf",
            results_path = config["paths"]["results"],
            f=["assembly_stats", "assembly_size_dist", "alignment_frequency"])

if config["assembly"]["metaspades"]:
    localrules:
        generate_metaspades_input
    rule generate_metaspades_input:
        """Generate input files for use with Metaspades"""
        input:
            lambda wildcards: get_all_assembly_files(assemblies[wildcards.assembly])
        output:
            R1=temp(expand("{results_path}/assembly/{{assembly}}/R1.fq",
                results_path=config["paths"]["results"])),
            R2=temp(expand("{results_path}/assembly/{{assembly}}/R2.fq",
                results_path=config["paths"]["results"])),
            se=touch(temp(expand("{results_path}/assembly/{{assembly}}/se.fq",
                results_path=config["paths"]["results"]))),
        params:
            assembly = lambda wildcards: assemblies[wildcards.assembly],
            assembler = "metaspades"
        script:
            "../scripts/assembly_utils.py"

    rule metaspades:
        input:
            R1=expand("{results_path}/assembly/{{assembly}}/R1.fq",
                results_path=config["paths"]["results"]),
            R2=expand("{results_path}/assembly/{{assembly}}/R2.fq",
                results_path=config["paths"]["results"]),
            se=expand("{results_path}/assembly/{{assembly}}/se.fq",
                results_path=config["paths"]["results"]),
        output:
            expand("{results_path}/assembly/{{assembly}}/final_contigs.fa",
                results_path=config["paths"]["results"])
        log:
            expand("{results_path}/assembly/{{assembly}}/spades.log",
                results_path=config["paths"]["results"])
        params:
            intermediate_contigs=expand("{results_path}/intermediate/assembly/{{assembly}}/intermediate_contigs",
                results_path=config["paths"]["results"]),
            corrected=expand("{results_path}/intermediate/assembly/{{assembly}}/corrected",
                results_path=config["paths"]["results"]),
            additional_settings=config["metaspades"]["extra_settings"],
            tmp=expand("{results_path}/{{assembly}}.metaspades",
                results_path=config["paths"]["temp"]),
            output_dir=expand("{results_path}/assembly/{{assembly}}",
                results_path=config["paths"]["results"])
        threads: config["metaspades"]["threads"]
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60*4
        conda:
            "../envs/metaspades.yml"
        shell:
            """
            # Create directories
            mkdir -p {params.tmp}
            # Clean output dir
            #rm -rf {params.output_dir}/*
            # Clean temp dir
            rm -rf {params.tmp}/*
            # Only use single-end if present
            if [ -s {input.se} ]; then
                single="-s {input.se}"
            else
                single=""
            fi
            metaspades.py \
                -t {threads} -1 {input.R1} -2 {input.R2} $single \
                -o {params.tmp} > {log} 2>&1

            # If set to keep intermediate contigs, move to intermediate folder before deleting
            if [ "{config[metaspades][keep_intermediate]}" == "True" ]; then
                mkdir -p {params.intermediate_contigs}
                cp -r {params.tmp}/K* {params.intermediate_contigs}
            fi
            if [ "{config[metaspades][keep_corrected]}" == "True" ]; then
                mkdir -p {params.corrected}
                cp -r {params.tmp}/corrected {params.corrected}
            fi

            # Clear intermediate contigs
            rm -rf {params.tmp}/K*
            # Clear corrected reads dir
            rm -rf {params.tmp}/corrected
            # Sync tmp output to outdir before removing
            cp -r {params.tmp}/* {params.output_dir}
            rm -rf {params.tmp}
            mv {params.output_dir}/scaffolds.fasta {params.output_dir}/final_contigs.fa
            """

else:
    localrules:
        generate_megahit_input
    rule generate_megahit_input:
        """Generate input lists for Megahit"""
        input:
            lambda wildcards: get_all_assembly_files(assemblies[wildcards.assembly])
        output:
            R1=temp(expand("{results_path}/assembly/{{assembly}}/input_1",
                results_path=config["paths"]["results"])),
            R2=temp(expand("{results_path}/assembly/{{assembly}}/input_2",
                results_path=config["paths"]["results"])),
            se=temp(expand("{results_path}/assembly/{{assembly}}/input_se",
                results_path=config["paths"]["results"]))
        log:
            expand("{results_path}/assembly/{{assembly}}/input_list.log",
                results_path=config["paths"]["results"])
        params:
            assembly = lambda wildcards: assemblies[wildcards.assembly]
        script:
            "../scripts/assembly_utils.py"

    rule megahit:
        input:
            R1=expand("{results_path}/assembly/{{assembly}}/input_1",
                results_path=config["paths"]["results"]),
            R2=expand("{results_path}/assembly/{{assembly}}/input_2",
                results_path=config["paths"]["results"]),
            se=expand("{results_path}/assembly/{{assembly}}/input_se",
                results_path=config["paths"]["results"]),
        output:
            expand("{results_path}/assembly/{{assembly}}/final_contigs.fa",
                results_path=config["paths"]["results"])
        log:
            expand("{results_path}/assembly/{{assembly}}/log",
                results_path=config["paths"]["results"])
        params:
            intermediate_contigs=expand("{results_path}/intermediate/assembly/{{assembly}}/intermediate_contigs",
                results_path=config["paths"]["results"]),
            additional_settings=config["megahit"]["extra_settings"],
            tmp=expand("{results_path}/{{assembly}}.megahit",
                results_path=config["paths"]["results"]),
            output_dir=lambda wildcards, output: os.path.dirname(output[0])
        threads: config["megahit"]["threads"]
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60*4
        conda:
            "../envs/megahit.yml"
        shell:
            """
            mkdir -p {config[paths][temp]}
            rm -rf {params.tmp}
            # Only use paired-end if present
            if [ -s {input.R1} ]; then
                R1=$(cat {input.R1})
                R2=$(cat {input.R2})
                paired="-1 $R1 -2 $R2"
            else
                paired=""
            fi
            # Only use single-end if present
            se=$(cat {input.se})
            if [ -s {input.se} ]; then
                single="-r $se"
            else
                single=""
            fi

            # Run Megahit
            megahit -t {threads} $paired $single -o {params.tmp} \
                {params.additional_settings} >{log} 2>&1

            # Sync intermediate contigs if asked for
            if [ "{config[megahit][keep_intermediate]}" == "True" ]; then
                mkdir -p {params.intermediate_contigs}
                cp -r {params.tmp}/intermediate_contigs/* {params.intermediate_contigs}
            fi

            # Cleanup intermediate
            rm -rf {params.tmp}/intermediate_contigs

            # Sync tmp output to outdir before removing
            cp -r {params.tmp}/* {params.output_dir}
            rm -rf {params.tmp}
            mv {params.output_dir}/final.contigs.fa {params.output_dir}/final_contigs.fa
            """

rule fasta2bed:
    """Creates bed-format file from assembly"""
    input:
        expand("{results_path}/assembly/{{assembly}}/final_contigs.fa",
            results_path=config["paths"]["results"])
    output:
        expand("{results_path}/assembly/{{assembly}}/final_contigs.bed",
            results_path=config["paths"]["results"])
    script:
        "../scripts/assembly_utils.py"

###########
# Mapping #
###########

rule bowtie_build:
    input:
        opj(config["paths"]["results"],"assembly","{assembly}","final_contigs.fa")
    output:
        expand(opj(config["paths"]["results"],"assembly","{{assembly}}","final_contigs.fa.{index}.bt2l"),index=range(1,5))
    params: prefix=opj(config["paths"]["results"],"assembly","{assembly}","final_contigs.fa")
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/quantify.yml"
    shell:
        """
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.prefix} \
            {params.prefix} > /dev/null 2>&1
        """

rule bowtie_map_pe:
    input:
        bt_index=expand(opj(config["paths"]["results"],"assembly","{{assembly}}",
                              "final_contigs.fa.{index}.bt2l"),
                          index=range(1,5)),
        R1=expand(opj(config["paths"]["results"], "intermediate", "preprocess",
                        "{{sample}}_{{unit}}_R1{p}.fastq.gz"), p=PREPROCESS),
        R2=expand(opj(config["paths"]["results"], "intermediate", "preprocess",
                        "{{sample}}_{{unit}}_R2{p}.fastq.gz"), p=PREPROCESS)
    output:
        bam=temp(opj(config["paths"]["results"],"assembly","{assembly}",
                     "mapping","{sample}_{unit}_pe.bam")),
        bai=temp(opj(config["paths"]["results"],"assembly","{assembly}",
                     "mapping","{sample}_{unit}_pe.bam.bai")),
        log=opj(config["paths"]["results"],"assembly","{assembly}",
                "mapping","{sample}_{unit}_pe.bam.log")
    params:
        temp_bam=opj(config["paths"]["temp"],"{assembly}-mapping-{sample}_{unit}_pe.bam"),
        setting=config["bowtie2"]["extra_settings"],
        prefix=opj(config["paths"]["results"],"assembly","{assembly}",
                     "final_contigs.fa")
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/quantify.yml"
    shell:
        """
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -1 {input.R1} -2 {input.R2} 2> {output.log} \
            | samtools view -bh - | samtools sort - -o {params.temp_bam}
        samtools index {params.temp_bam}
        mv {params.temp_bam} {output.bam}
        mv {params.temp_bam}.bai {output.bai}
        """

rule bowtie_map_se:
    input:
        bt_index=expand(opj(config["paths"]["results"],"assembly","{{assembly}}",
                              "final_contigs.fa.{index}.bt2l"),
                          index=range(1,5)),
        se=expand(opj(config["paths"]["results"], "intermediate", "preprocess",
                        "{{sample}}_{{unit}}_se{p}.fastq.gz"), p=PREPROCESS)
    output:
        bam=temp(opj(config["paths"]["results"],"assembly","{assembly}",
                     "mapping","{sample}_{unit}_se.bam")),
        bai=temp(opj(config["paths"]["results"],"assembly","{assembly}",
                     "mapping","{sample}_{unit}_se.bam.bai")),
        log=opj(config["paths"]["results"],"assembly","{assembly}",
                "mapping","{sample}_{unit}_se.bam.log")
    params:
        temp_bam=opj(config["paths"]["temp"],"{assembly}-mapping-{sample}_{unit}_se.bam"),
        setting=config["bowtie2"]["extra_settings"],
        prefix=opj(config["paths"]["results"],"assembly","{assembly}",
                     "final_contigs.fa")
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/quantify.yml"
    shell:
        """
        bowtie2 {params.setting} -p {threads} -x {params.prefix} \
            -U {input.se} 2>{output.log} | samtools view -bh - | samtools sort - -o {params.temp_bam}
        samtools index {params.temp_bam}
        mv {params.temp_bam} {output.bam}
        mv {params.temp_bam}.bai {output.bai}
        """

##############
# Statistics #
##############

rule samtools_flagstat:
    """
    Generate mapping statistics
    """
    input:
        lambda wildcards: get_bamfiles(wildcards.assembly,
                                       assemblies[wildcards.assembly],
                                       config["paths"]["results"], POSTPROCESS)
    output:
        opj(config["paths"]["results"],"assembly","{assembly}",
                 "mapping","flagstat.tsv")
    params:
        post = POSTPROCESS
    conda:
        "../envs/quantify.yml"
    shell:
        """
        for f in {input} ;
        do
            al=$(samtools \
                flagstat \
                $f | grep " mapped (" | cut -f2 -d '('| cut -f1 -d ' ')
            n=$(basename $f | sed 's/_[ps]e{params.post}.bam//g')
            echo -e "$n\t$al" >> {output}
        done
        """

rule assembly_stats:
    input:
        fa = expand(opj(config["paths"]["results"],"assembly","{assembly}",
                   "final_contigs.fa"), assembly=assemblies.keys()),
        flagstat = expand(opj(config["paths"]["results"],"assembly","{assembly}",
                   "mapping","flagstat.tsv"), assembly = assemblies.keys())
    output:
        report(opj(config["paths"]["results"], "report", "assembly", "assembly_stats.tsv"),
               caption="../report/assembly.rst", category="Assembly"),
        opj(config["paths"]["results"], "report", "assembly", "assembly_size_dist.tsv")
    script:
        "../scripts/assembly_utils.py"

rule plot_assembly_stats:
    input:
        stat = opj(config["paths"]["results"], "report", "assembly", "assembly_stats.tsv"),
        dist = opj(config["paths"]["results"], "report", "assembly", "assembly_size_dist.tsv"),
        maps = expand(opj(config["paths"]["results"],"assembly","{assembly}",
                 "mapping","flagstat.tsv"), assembly = assemblies.keys())
    output:
        report(opj(config["paths"]["results"], "report", "assembly",
            "assembly_stats.pdf"), caption="../report/assembly.rst",
               category="Assembly"),
        report(opj(config["paths"]["results"], "report", "assembly",
                   "assembly_size_dist.pdf"), caption="../report/assembly.rst",
               category="Assembly"),
        report(opj(config["paths"]["results"], "report", "assembly",
                   "alignment_frequency.pdf"), caption="../report/assembly.rst",
               category="Assembly")
    conda:
        "../envs/plotting.yml"
    notebook:
        "../notebooks/assembly_stats.py.ipynb"
