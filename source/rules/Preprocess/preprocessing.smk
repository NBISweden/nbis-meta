localrules:
    link_files,
    sortmerna_unzip_fastq,
    sortmerna_zip_aligned_fastq,
    sortmerna_zip_other_fastq,
    sortmerna_link_pe,
    sortmerna_link_se,
    download_phix,
    bowtie_build_phix,
    fastuniq_se,
    avg_seq_length

def link(target,link_name):
    target_abs = os.path.abspath(target)
    link_abs = os.path.abspath(link_name)
    shell("ln -s {target_abs} {link_abs}")

def get_interleaved(sample,runID):
    files = []
    if "interleaved" in samples[sample][runID].keys():
        inter = samples[sample][runID]["interleaved"]
        R1 = samples[sample][runID]["R1"]
        R2 = samples[sample][runID]["R2"]
        files.append(inter)
    else:
        files.append("")
    return files

rule deinterleave_fastq:
    input:
        lambda wildcards: get_interleaved(wildcards.sample,wildcards.run)
    output:
        R1=opj(config["intermediate_path"],"deinterleaved","{sample}_{run}_R1.fastq.gz"),
        R2=opj(config["intermediate_path"],"deinterleaved","{sample}_{run}_R2.fastq.gz")
    params:
        script="source/utils/deinterleave_fastq.sh",
        tmp_r1=opj(os.path.expandvars(config["scratch_path"]),"{sample}_{run}_R1.fastq.gz"),
        tmp_r2=opj(os.path.expandvars(config["scratch_path"]),"{sample}_{run}_R2.fastq.gz")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*24
    run:
        for item in input:
            if not item:
                continue
            shell("echo {params.script} {item} {params.tmp_r1} {params.tmp_r2} compress")
            shell("{params.script} {item} {params.tmp_r1} {params.tmp_r2} compress")
            shell("mv {params.tmp_r1} {output.R1}")
            shell("mv {params.tmp_r2} {output.R2}")

rule link_files:
    """Symlink sample to make downstream processing easier"""
    input:
        lambda wildcards: samples[wildcards.sample][wildcards.run][wildcards.pair]
    output:
        opj(config["intermediate_path"],"preprocess","{sample}_{run}_{pair}.fastq.gz")
    message: "Linking {wildcards.sample}_{wildcards.run}_{wildcards.pair}.fastq.gz"
    run:
        cmd="ln -s {} {}".format(
            os.path.abspath(
                samples[wildcards.sample][wildcards.run][wildcards.pair]),
                str(output))
        shell(cmd)

rule sortmerna_merge_fastq:
    """Merge fastq output from SortMeRNA"""
    input:
        R1=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R1.fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R2.fastq.gz")
    output:
        temp(opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_merged.fastq"))
    params:
        scratch=os.path.expandvars(config["scratch_path"]),
        R1_unzipped=opj(os.path.expandvars(config["scratch_path"]),
                          "{sample}_{run}_R1.fastq"),
        R2_unzipped=opj(os.path.expandvars(config["scratch_path"]),
                        "{sample}_{run}_R2.fastq"),
        merged=opj(os.path.expandvars(config["scratch_path"]),
                   "{sample}_{run}_merged.fastq")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        mkdir -p {params.scratch}
        # Unzip to scratch dir
        gunzip -c {input.R1} > {params.R1_unzipped}
        gunzip -c {input.R2} > {params.R2_unzipped}
        # Merge
        merge-paired-reads.sh \
            {params.R1_unzipped} \
            {params.R2_unzipped} \
            {params.merged} >/dev/null 2>&1
        # Move output
        mv {params.merged} {output[0]}
        # Clean up
        rm {params.R1_unzipped} {params.R2_unzipped}
        """

def get_sortmerna_ref_string(path, s):
    """
    Constructs the SortMeRNA --ref string

    :param path: Resource path from config
    :param s: Sortmerna databases from config
    :return: STRING,STRING formatted string
    """
    files=["{p}/rRNA_databases/{db}".format(p=path, db=db)
           for db in config["sortmerna_dbs"]]
    ref_string =":".join(
        ["{},{}".format(f,f) for f in files])
    return ref_string

rule sortmerna_fastq_pe:
    """Run SortMeRNA on paired end input"""
    input:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_merged.fastq"),
        expand(opj(config["resource_path"],"rRNA_databases","{file}.{suffix}"),
            suffix=["bursttrie_0.dat","kmer_0.dat","pos_0.dat","stats"],
            file=config["sortmerna_dbs"])
    output:
        aligned=temp(opj(config["intermediate_path"],"preprocess",
                         "{sample}_{run}_merged.rRNA.fastq")),
        other=temp(opj(config["intermediate_path"],"preprocess",
                       "{sample}_{run}_merged.non_rRNA.fastq"))
    log:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_pe.sortmerna.log")
    params:
        paired_strategy=config["sortmerna_paired_strategy"],
        score_params=config["sortmerna_params"],
        other_prefix=opj(config["scratch_path"],"{sample}_{run}_merged.non_rRNA"),
        aligned_prefix=opj(config["scratch_path"],"{sample}_{run}_merged.rRNA"),
        scratch=config["scratch_path"],
        ref_string=get_sortmerna_ref_string(config["resource_path"],
                                            config["sortmerna_dbs"])
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    shell:
         """
         mkdir -p {params.scratch} 
         # Run SortMeRNA
         sortmerna \
            --blast 1 \
            --log \
            -v \
            --fastx \
            --ref {params.ref_string} \
            --reads {input[0]} \
            -a {threads} \
            --{params.paired_strategy} \
            --aligned {params.aligned_prefix} \
            --other {params.other_prefix} \
            {params.score_params} \
            >/dev/null 2>&1
         
         mv {params.aligned_prefix}.fastq {output.aligned}
         mv {params.aligned_prefix}.log {log}
         mv {params.other_prefix}.fastq {output.other}
         """

rule sortmerna_split_rRNA_fastq:
    input:
        aligned=opj(config["intermediate_path"],"preprocess",
                      "{sample}_{run}_merged.rRNA.fastq"),
    output:
        R1=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R1.rRNA.fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_R2.rRNA.fastq.gz")
    params:
        tmpdir=opj(os.path.expandvars(config["scratch_path"]),
                   "{sample}_{run}_sortmerna"),
        R1=opj(os.path.expandvars(config["scratch_path"]),
                 "{sample}_{run}_sortmerna","{sample}_{run}_R1.rRNA.fastq"),
        R2=opj(os.path.expandvars(config["scratch_path"]),
                 "{sample}_{run}_sortmerna","{sample}_{run}_R2.rRNA.fastq")
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        unmerge-paired-reads.sh \
            {input.aligned} \
            {params.R1} \
            {params.R2} >/dev/null 2>&1
        gzip {params.R1}
        gzip {params.R2}
        mv {params.R1}.gz {output.R1}
        mv {params.R2}.gz {output.R2}
        """

rule sortmerna_split_other_fastq:
    input:
        other=opj(config["intermediate_path"],"preprocess",
                  "{sample}_{run}_merged.non_rRNA.fastq")
    output:
        R1=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R1.non_rRNA.fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
               "{sample}_{run}_R2.non_rRNA.fastq.gz")
    params:
        tmpdir=opj(os.path.expandvars(config["scratch_path"]),
                   "{sample}_{run}_sortmerna"),
        R1=opj(os.path.expandvars(config["scratch_path"]),
               "{sample}_{run}_sortmerna","{sample}_{run}_R1.non_rRNA.fastq"),
        R2=opj(os.path.expandvars(config["scratch_path"]),
               "{sample}_{run}_sortmerna","{sample}_{run}_R2.non_rRNA.fastq")
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        mkdir -p {params.tmpdir}
        unmerge-paired-reads.sh \
            {input.other} \
            {params.R1} \
            {params.R2} >/dev/null 2>&1
        gzip {params.R1}
        gzip {params.R2}
        mv {params.R1}.gz {output.R1}
        mv {params.R2}.gz {output.R2}
        """

rule sortmerna_unzip_fastq:
    input:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se.fastq.gz")
    output:
        temp(opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_se.fastq"))
    shell:
        """
        gunzip -c {input[0]} > {output[0]}
        """

rule sortmerna_fastq_se:
    input:
        opj(config["intermediate_path"],"preprocess",
                  "{sample}_{run}_se.fastq"),
        expand(opj(config["resource_path"],"rRNA_databases","{file}.{suffix}"),
            suffix=["bursttrie_0.dat","kmer_0.dat","pos_0.dat","stats"],
            file=config["sortmerna_dbs"])
    output:
        aligned=temp(opj(config["intermediate_path"],"preprocess",
                         "{sample}_{run}_se.rRNA.fastq")),
        other=temp(opj(config["intermediate_path"],"preprocess",
                       "{sample}_{run}_se.non_rRNA.fastq"))
    log:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se.sortmerna.log")
    params:
        score_params=config["sortmerna_params"],
        other_prefix=opj(config["scratch_path"],"{sample}_{run}_se.non_rRNA"),
        aligned_prefix=opj(config["scratch_path"],"{sample}_{run}_se.rRNA"),
        scratch=config["scratch_path"],
        ref_string=get_sortmerna_ref_string(config["resource_path"],
                                            config["sortmerna_dbs"])
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        mkdir -p {params.scratch}
        # Run SortMeRNA
        sortmerna \
            --blast 1 \
            --log \
            -v \
            --fastx \
            --ref {params.ref_string} \
            --reads {input[0]} \
            -a {threads} \
            --aligned {params.aligned_prefix} \
            --other {params.other_prefix} \
            {params.score_params} \
            >/dev/null 2>&1
        
        mv {params.aligned_prefix}.fastq {output.aligned}
        mv {params.aligned_prefix}.log {log}
        mv {params.other_prefix}.fastq {output.other}
        """

rule sortmerna_zip_aligned_fastq:
    input:
        fastq=opj(config["intermediate_path"],"preprocess",
                  "{sample}_{run}_se.rRNA.fastq")
    output:
        fastq=opj(config["intermediate_path"],"preprocess",
                    "{sample}_{run}_se.rRNA.fastq.gz")
    shell:
        """
        gzip {input.fastq}
        """

rule sortmerna_zip_other_fastq:
    input:
        fastq=opj(config["intermediate_path"],"preprocess",
                  "{sample}_{run}_se.non_rRNA.fastq")
    output:
        fastq=opj(config["intermediate_path"],"preprocess",
                    "{sample}_{run}_se.non_rRNA.fastq.gz")
    shell:
        """
        gzip {input.fastq}
        """

rule sortmerna_link_pe:
  input:
    R1=opj(config["intermediate_path"],"preprocess",
                "{sample}_{run}_R1."+config["sortmerna_keep"]+".fastq.gz"),
    R2=opj(config["intermediate_path"],"preprocess",
             "{sample}_{run}_R2."+config["sortmerna_keep"]+".fastq.gz")
  output:
    R1=opj(config["intermediate_path"],"preprocess",
             "{sample}_{run}_R1.sortmerna.fastq.gz"),
    R2=opj(config["intermediate_path"],"preprocess",
           "{sample}_{run}_R2.sortmerna.fastq.gz")
  run:
    link(input.R1, output.R1)
    link(input.R2, output.R2)

rule sortmerna_link_se:
    input:
        se=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_se."+config["sortmerna_keep"]+".fastq.gz")
    output:
        se=opj(config["intermediate_path"],"preprocess",
                 "{sample}_{run}_se.sortmerna.fastq.gz")
    run:
        link(input.se, output.se)

def get_trimmomatic_string(seq_type):
    """
    Generates trimsetting string for Trimmomatic

    :param seq_type: PE or SE depending on sequencing type
    :return: trimsettings string
    """
    trim_adapters=config["trim_adapters"]
    adapter_fasta_dir="$CONDA_PREFIX/share/trimmomatic/adapters"
    adapter="{}/{}.fa".format(adapter_fasta_dir,
                            config["trimmomatic_{}_adapter".format(seq_type)])
    adapter_params=config["{}_adapter_params".format(seq_type)]
    pre_adapter_params=config["{}_pre_adapter_params".format(seq_type)]
    post_adapter_params=config["{}_post_adapter_params".format(seq_type)]
    trimsettings=pre_adapter_params
    if trim_adapters:
        trimsettings+=" ILLUMINACLIP:"+adapter+":"+adapter_params
    trimsettings+=" "+post_adapter_params
    return trimsettings

rule trimmomatic_pe:
    input:
        R1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["trimming"]+".fastq.gz")
    output:
        R1P=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".trimmomatic.fastq.gz"),
        R1U=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".trimmomatic.U.fastq.gz"),
        R2P=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["trimming"]+".trimmomatic.fastq.gz"),
        R2U=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["trimming"]+".trimmomatic.U.fastq.gz"),
        R1log=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".trimmomatic.log"),
        R2log=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["trimming"]+".trimmomatic.log")
    params:
        trim_string=get_trimmomatic_string("pe")
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input.R1} {input.R2} \
            {output.R1P} {output.R1U} \
            {output.R2P} {output.R2U} \
            {params.trim_string} \
            2>{output.R1log}
        sed \
            's/{wildcards.sample}_{wildcards.run}_R1/{wildcards.sample}_{wildcards.run}_R2/g' \
            {output.R1log} > {output.R2log}    
        """

rule trimmomatic_se:
    """Run Trimmomatic on single-end input"""
    input:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["trimming"]+".fastq.gz")
    output:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["trimming"]+".trimmomatic.fastq.gz"),
    log:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["trimming"]+".trimmomatic.log")
    params:
        trim_string=get_trimmomatic_string("se")
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        trimmomatic SE \
            -threads {threads} \
            {input} \
            {output} \
            {params.trim_string} \
            2>{log}
        """

rule cutadapt_pe:
    input:
        R1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["trimming"]+".fastq.gz")
    output:
        fastq1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".cutadapt.fastq.gz"),
        fastq2=opj(config["intermediate_path"],
            "preprocess","{sample}_{run}_R2"+preprocess_suffices["trimming"]+".cutadapt.fastq.gz"),
        R1log=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["trimming"]+".cutadapt.log"),
        R2log=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["trimming"]+".cutadapt.log")
    params:
        adapter=config["adapter_sequence"],
        rev_adapter=config["rev_adapter_sequence"],
        error_rate=config["cutadapt_error_rate"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    threads: 10
    shell:
        """
        cutadapt \
            -e {params.error_rate} \
            -a {params.adapter} \
            -A {params.rev_adapter} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            -j {threads} \
            {input.R1} {input.R2} > {output.R1log}
        cp {output.R1log} {output.R2log}
        """

rule cutadapt_se:
    input:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["trimming"]+".fastq.gz")
    output:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["trimming"]+".cutadapt.fastq.gz")
    log:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["trimming"]+".cutadapt.log")
    params:
        adapter=config["adapter_sequence"],
        error_rate=config["cutadapt_error_rate"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    threads: 10
    shell:
        """
        cutadapt \
            -e {params.error_rate} \
            -a {params.adapter} \
            -j {threads} \
            -o {output[0]} {input[0]} > {log}
        """

rule download_phix:
    """Downloads the phiX genome"""
    output:
        opj(config["resource_path"],"phix","phix.fasta")
    params:
        url_base="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015"
    shell:
        """
        curl \
            -L \
            -o {output[0]}.gz \
            -s \
            {params.url_base}/GCF_000819615.1_ViralProj14015_genomic.fna.gz
        gunzip {output[0]}.gz
        """

rule bowtie_build_phix:
    """Build bowtie2 index for phiX"""
    input:
        opj(config["resource_path"],"phix","phix.fasta")
    output:
        expand(opj(config["resource_path"],"phix","phix.{index}.bt2"),
               index=range(1,5))
    params:
        prefix=opj(config["resource_path"],"phix","phix")
    threads: 1
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input} {params.prefix} >/dev/null 2>&1
        """

rule filter_phix_pe:
    """Maps reads against the phiX genome, keeping non-concordantly mapped"""
    input:
        bt_index=expand(opj(config["resource_path"],"phix",
                            "phix.{index}.bt2"),index=range(1,5)),
        R1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["phixfilt"]+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["phixfilt"]+".fastq.gz")
    output:
        R1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["phixfilt"]+".phixfilt.fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["phixfilt"]+".phixfilt.fastq.gz"),
        log=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_PHIX_pe"+preprocess_suffices["phixfilt"]+".log")
    params:
        tmp_out=config["scratch_path"],
        setting=config["bowtie2_params"],
        prefix=opj(config["resource_path"],"phix","phix")
    threads: config["bowtie2_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        mkdir -p {params.tmp_out}
        bowtie2 \
            {params.setting} \
            -p {threads} \
            -x {params.prefix} \
            -1 {input.R1} \
            -2 {input.R2} \
            --un-conc-gz \
            {params.tmp_out}/{wildcards.sample}_{wildcards.run}_R%.filtered.fastq.gz > /dev/null 2>{output.log}
        mv {params.tmp_out}/{wildcards.sample}_{wildcards.run}_R1.filtered.fastq.gz {output.R1}
        mv {params.tmp_out}/{wildcards.sample}_{wildcards.run}_R2.filtered.fastq.gz {output.R2}
        """

rule filter_phix_se:
    """Maps reads against the phiX genome, keeping non-concordantly mapped"""
    input:
        bt_index=expand(opj(config["resource_path"],"phix",
                            "phix.{index}.bt2"),index=range(1,5)),
        se=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["phixfilt"]+".fastq.gz")
    output:
        se=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["phixfilt"]+".phixfilt.fastq.gz"),
    log:
        opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_PHIX_se"+preprocess_suffices["phixfilt"]+".log")
    params:
        tmp_out=config["scratch_path"],
        setting=config["bowtie2_params"],
        prefix=opj(config["resource_path"],"phix","phix")
    threads: config["bowtie2_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        mkdir -p {params.tmp_out}
        bowtie2 \
            {params.setting} \
            -p {threads} \
            -x {params.prefix} \
            --un-gz \
            {params.tmp_out}/{wildcards.sample}_{wildcards.run}_se.filtered.fastq.gz {input.se} > /dev/null 2>{log}
        mv {params.tmp_out}/{wildcards.sample}_{wildcards.run}_se.filtered.fastq.gz {output.se}
        """

rule fastuniq:
    input:
        R1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["fastuniq"]+".fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["fastuniq"]+".fastq.gz")
    output:
        R1=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R1"+preprocess_suffices["fastuniq"]+".fastuniq.fastq.gz"),
        R2=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_R2"+preprocess_suffices["fastuniq"]+".fastuniq.fastq.gz"),
    params:
        R1_intmp=opj(config["scratch_path"],
            "{sample}_{run}_R1"+preprocess_suffices["fastuniq"]+".fastq"),
        R2_intmp=opj(config["scratch_path"],
            "{sample}_{run}_R2"+preprocess_suffices["fastuniq"]+".fastq"),
        R1_outtmp=opj(config["scratch_path"],
            "{sample}_{run}_R1"+preprocess_suffices["fastuniq"]+".fastuniq.fastq"),
        R2_outtmp=opj(config["scratch_path"],
            "{sample}_{run}_R2"+preprocess_suffices["fastuniq"]+".fastuniq.fastq"),
        file_list=opj(config["scratch_path"],
            "{sample}_{run}.filelist")
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../../../envs/preprocess.yml"
    shell:
        """
        gunzip -c {input.R1} > {params.R1_intmp}
        gunzip -c {input.R2} > {params.R2_intmp}
        echo {params.R1_intmp} > {params.file_list}
        echo {params.R2_intmp} >> {params.file_list}
        fastuniq \
            -i {params.file_list} \
            -t q \
            -o {params.R1_outtmp} \
            -p {params.R2_outtmp}
        gzip -c {params.R1_outtmp} > {output.R1}
        gzip -c {params.R2_outtmp} > {output.R2}
        """

rule fastuniq_se:
    """Dummy rule for fastuniq on single-end input"""
    input:
        se=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["fastuniq"]+".fastq.gz")
    output:
        se=opj(config["intermediate_path"],"preprocess",
            "{sample}_{run}_se"+preprocess_suffices["fastuniq"]+".fastuniq.fastq.gz")
    shell:
        """
        mv {input.se} {output.se}
        """

rule avg_seq_length:
    """Extracts average sequence lengths from FastQC"""
    input:
        "results/report/samples_report_data/multiqc_general_stats.txt"
    output:
        opj(config["intermediate_path"],"preprocess","read_lengths.tsv")
    run:
        import pandas as pd
        df=pd.read_csv(input[0], header=0, sep="\t", index_col=0)
        df=df.loc[:,"FastQC_mqc-generalstats-fastqc-avg_sequence_length"]
        df.columns=["read_length"]
        df.to_csv(output[0], sep="\t", index=True)