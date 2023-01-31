localrules:
    link_files,
    download_rRNA_database,
    sortmerna_unzip_fastq,
    sortmerna_zip_aligned_fastq,
    sortmerna_zip_other_fastq,
    sortmerna_link_pe,
    sortmerna_link_se,
    download_phix,
    bowtie_build_phix,
    fastuniq_se,
    multiqc,


##### master rule for preprocessing #####


rule qc:
    input:
        preprocessing_input(config),


##### utility rules #####


rule link_files:
    """Symlink sample to make downstream processing easier"""
    input:
        lambda wildcards: samples[wildcards.sample][wildcards.unit][wildcards.pair],
    output:
        results + "/intermediate/preprocess/{sample}_{unit}_{pair}.fastq.gz",
    run:
        link(input[0], output[0])


##### sortmerna #####


rule download_rRNA_database:
    output:
        "resources/rRNA_databases/{file}.fasta",
    log:
        "resources/rRNA_databases/{file}.dl.log",
    retries: 3
    params:
        url="https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/{file}.fasta",
    shell:
        """
         curl -L -o {output} {params.url} > {log} 2>&1
         """


rule index_db:
    input:
        fasta="resources/rRNA_databases/{file}.fasta",
    output:
        expand(
            "resources/rRNA_databases/{{file}}.fasta.{suffix}",
            suffix=["bursttrie_0.dat", "kmer_0.dat", "pos_0.dat", "stats"],
        ),
    log:
        "resources/rRNA_databases/{file}.index.log",
    resources:
        runtime=300,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "SortMeRNA/2.1b",
    shell:
        """
        indexdb_rna --ref {input.fasta},{input.fasta} > {log} 2>&1
        """


rule sortmerna_merge_fastq:
    """Merge fastq output from SortMeRNA"""
    input:
        R1=results + "/intermediate/preprocess/{sample}_{unit}_R1.fastq.gz",
        R2=results + "/intermediate/preprocess/{sample}_{unit}_R2.fastq.gz",
    output:
        temp(results + "/intermediate/preprocess/{sample}_{unit}_merged.fastq"),
    log:
        results + "/intermediate/preprocess/{sample}_{unit}.sortmerna_merge.log",
    params:
        scratch=temppath,
        R1_unzipped=results + "/{sample}_{unit}_R1.fastq",
        R2_unzipped=results + "/{sample}_{unit}_R2.fastq",
        merged=temppath + "/{sample}_{unit}_merged.fastq",
    resources:
        runtime=360,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    group:
        "sortmerna"
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "SortMeRNA/2.1b",
    shell:
        """
        mkdir -p {params.scratch}
        # Unzip to scratch dir
        gunzip -c {input.R1} > {params.R1_unzipped}
        gunzip -c {input.R2} > {params.R2_unzipped}
        # Merge
        merge-paired-reads.sh {params.R1_unzipped} {params.R2_unzipped} \
            {params.merged} >{log} 2>&1
        # Move output
        mv {params.merged} {output}
        # Clean up
        rm {params.R1_unzipped} {params.R2_unzipped}
        """


rule sortmerna_fastq_pe:
    """Run SortMeRNA on paired end input"""
    input:
        fastq=results + "/intermediate/preprocess/{sample}_{unit}_merged.fastq",
        db=expand(
            "resources/rRNA_databases/{file}.{suffix}",
            suffix=["bursttrie_0.dat", "kmer_0.dat", "pos_0.dat", "stats"],
            file=config["sortmerna"]["dbs"],
        ),
    output:
        aligned=temp(
            results + "/intermediate/preprocess/{sample}_{unit}_merged.rRNA.fastq"
        ),
        other=temp(
            results + "/intermediate/preprocess/{sample}_{unit}_merged.non_rRNA.fastq"
        ),
    log:
        results + "/intermediate/preprocess/{sample}_{unit}_pe.sortmerna.log",
    params:
        paired_strategy=config["sortmerna"]["paired_strategy"],
        score_params=config["sortmerna"]["extra_settings"],
        other_prefix=temppath + "/{sample}_{unit}_merged.non_rRNA",
        aligned_prefix=temppath + "/{sample}_{unit}_merged.rRNA",
        scratch=temppath,
        ref_string=get_sortmerna_ref_string(config["sortmerna"]["dbs"]),
    threads: 10
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    group:
        "sortmerna"
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "SortMeRNA/2.1b",
    shell:
        """
         mkdir -p {params.scratch}
         # Run SortMeRNA
         sortmerna --blast 1 --log -v --fastx --ref {params.ref_string} \
            --reads {input.fastq} -a {threads} --{params.paired_strategy} \
            --aligned {params.aligned_prefix} --other {params.other_prefix} \
            {params.score_params} >{log} 2>&1

         mv {params.aligned_prefix}.fastq {output.aligned}
         mv {params.aligned_prefix}.log {log}
         mv {params.other_prefix}.fastq {output.other}
         """


rule sortmerna_split_fastq:
    input:
        other=results
        + "/intermediate/preprocess/{sample}_{unit}_merged.{RNA_type}.fastq",
    output:
        R1=results + "/intermediate/preprocess/{sample}_{unit}_R1.{RNA_type}.fastq.gz",
        R2=results + "/intermediate/preprocess/{sample}_{unit}_R2.{RNA_type}.fastq.gz",
    log:
        results
        + "/intermediate/preprocess/{sample}_{unit}.sortmerna_unmerge.{RNA_type}.log",
    params:
        tmpdir=temppath + "/{sample}_{unit}_sortmerna",
        R1=temppath + "/{sample}_{unit}_sortmerna/{sample}_{unit}_R1.{RNA_type}.fastq",
        R2=temppath + "/{sample}_{unit}_sortmerna/{sample}_{unit}_R2.{RNA_type}.fastq",
    resources:
        runtime=360,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    group:
        "sortmerna"
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "SortMeRNA/2.1b",
    shell:
        """
        mkdir -p {params.tmpdir}
        unmerge-paired-reads.sh {input.other} {params.R1} {params.R2} >{log} 2>&1
        gzip {params.R1}
        gzip {params.R2}
        mv {params.R1}.gz {output.R1}
        mv {params.R2}.gz {output.R2}
        """


rule sortmerna_unzip_fastq:
    input:
        results + "/intermediate/preprocess/{sample}_{unit}_se.fastq.gz",
    output:
        temp(results + "/intermediate/preprocess/{sample}_{unit}_se.fastq"),
    log:
        results + "/intermediate/preprocess/{sample}_{unit}_se.sortmerna_unzip.log",
    shell:
        """
        gunzip -c {input} > {output} 2>{log}
        """


rule sortmerna_fastq_se:
    input:
        fastq=results + "/intermediate/preprocess/{sample}_{unit}_se.fastq",
        db=expand(
            "resources/rRNA_databases/{file}.{suffix}",
            suffix=["bursttrie_0.dat", "kmer_0.dat", "pos_0.dat", "stats"],
            file=config["sortmerna"]["dbs"],
        ),
    output:
        aligned=temp(results + "/intermediate/preprocess/{sample}_{unit}_se.rRNA.fastq"),
        other=temp(
            results + "/intermediate/preprocess/{sample}_{unit}_se.non_rRNA.fastq"
        ),
    log:
        results + "/intermediate/preprocess/{sample}_{unit}_se.sortmerna.log",
    params:
        score_params=config["sortmerna"]["extra_settings"],
        other_prefix=temppath + "/{sample}_{unit}_se.non_rRNA",
        aligned_prefix=temppath + "/{sample}_{unit}_se.rRNA",
        scratch=temppath,
        ref_string=get_sortmerna_ref_string(config["sortmerna"]["dbs"]),
    threads: 10
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "SortMeRNA/2.1b",
    shell:
        """
        mkdir -p {params.scratch}
        # Run SortMeRNA
        sortmerna --blast 1 --log -v --fastx --ref {params.ref_string} \
            --reads {input.fastq} -a {threads} --other {params.other_prefix} \
            --aligned {params.aligned_prefix} {params.score_params} \
            >{log} 2>&1

        mv {params.aligned_prefix}.fastq {output.aligned}
        mv {params.aligned_prefix}.log {log}
        mv {params.other_prefix}.fastq {output.other}
        """


rule sortmerna_zip_aligned_fastq:
    input:
        fastq=results + "/intermediate/preprocess/{sample}_{unit}_se.rRNA.fastq",
    output:
        fastq=results + "/intermediate/preprocess/{sample}_{unit}_se.rRNA.fastq.gz",
    log:
        results + "/intermediate/preprocess/{sample}_{unit}_se.sortmerna_zip_rRNA.log",
    shell:
        """
        gzip {input.fastq} 2>{log}
        """


rule sortmerna_zip_other_fastq:
    input:
        fastq=results + "/intermediate/preprocess/{sample}_{unit}_se.non_rRNA.fastq",
    output:
        fastq=results + "/intermediate/preprocess/{sample}_{unit}_se.non_rRNA.fastq.gz",
    log:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se.sortmerna_zip_non_rRNA.log",
    shell:
        """
        gzip {input.fastq} 2>{log}
        """


rule sortmerna_link_pe:
    input:
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1."
        + config["sortmerna"]["keep"]
        + ".fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2."
        + config["sortmerna"]["keep"]
        + ".fastq.gz",
    output:
        R1=results + "/intermediate/preprocess/{sample}_{unit}_R1.sortmerna.fastq.gz",
        R2=results + "/intermediate/preprocess/{sample}_{unit}_R2.sortmerna.fastq.gz",
    run:
        link(input.R1, output.R1)
        link(input.R2, output.R2)


rule sortmerna_link_se:
    input:
        se=results
        + "/intermediate/preprocess/{sample}_{unit}_se."
        + config["sortmerna"]["keep"]
        + ".fastq.gz",
    output:
        se=results + "/intermediate/preprocess/{sample}_{unit}_se.sortmerna.fastq.gz",
    run:
        link(input.se, output.se)


##### trimmomatic #####


rule trimmomatic_pe:
    input:
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".fastq.gz",
    output:
        R1P=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.fastq.gz",
        R1U=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.U.fastq.gz",
        R2P=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.fastq.gz",
        R2U=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.U.fastq.gz",
    log:
        R1log=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.log",
        R2log=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.log",
    params:
        trim_string=get_trimmomatic_string("pe", config),
    threads: 10
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "trimmomatic/0.39",
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input.R1} {input.R2} \
            {output.R1P} {output.R1U} \
            {output.R2P} {output.R2U} \
            {params.trim_string} \
            2>{log.R1log}
        sed \
            's/{wildcards.sample}_{wildcards.unit}_R1/{wildcards.sample}_{wildcards.unit}_R2/g' \
            {log.R1log} > {log.R2log}
        """


rule trimmomatic_se:
    """Run Trimmomatic on single-end input"""
    input:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["trimming"]
        + ".fastq.gz",
    output:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.fastq.gz",
    log:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["trimming"]
        + ".trimmomatic.log",
    params:
        trim_string=get_trimmomatic_string("se", config),
    threads: 10
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "trimmomatic/0.39",
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
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".fastq.gz",
    output:
        fastq1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".cutadapt.fastq.gz",
        fastq2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".cutadapt.fastq.gz",
    log:
        R1log=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["trimming"]
        + ".cutadapt.log",
        R2log=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".cutadapt.log",
        err=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["trimming"]
        + ".cutadapt.err",
    params:
        adapter=config["cutadapt"]["adapter_sequence"],
        rev_adapter=config["cutadapt"]["rev_adapter_sequence"],
        error_rate=config["cutadapt"]["error_rate"],
        extra_params=config["cutadapt"]["extra_params"],
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "cutadapt/2.3",
    threads: 10
    shell:
        """
        cutadapt \
            {params.extra_params} \
            -e {params.error_rate} \
            -a {params.adapter} \
            -A {params.rev_adapter} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            -j {threads} \
            {input.R1} {input.R2} > {log.R1log} 2>{log.err}
        cp {log.R1log} {log.R2log}
        """


rule cutadapt_se:
    input:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["trimming"]
        + ".fastq.gz",
    output:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["trimming"]
        + ".cutadapt.fastq.gz",
    log:
        results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["trimming"]
        + ".cutadapt.log",
    params:
        adapter=config["cutadapt"]["adapter_sequence"],
        error_rate=config["cutadapt"]["error_rate"],
        extra_params=config["cutadapt"]["extra_params"],
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "cutadapt/2.3",
    threads: 10
    shell:
        """
        cutadapt \
            {params.extra_params} \
            -e {params.error_rate} \
            -a {params.adapter} \
            -j {threads} \
            -o {output} {input} > {log}
        """


rule download_phix:
    """Downloads the phiX genome"""
    output:
        "resources/phix/phix.fasta",
    log:
        "resources/phix/phix.log",
    retries: 3
    params:
        url_base="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015",
    shell:
        """
        curl \
            -L \
            -o {output}.gz \
            -s \
            {params.url_base}/GCF_000819615.1_ViralProj14015_genomic.fna.gz \
            > {log} 2>&1
        gunzip {output}.gz
        """


rule bowtie_build_phix:
    """Build bowtie2 index for phiX"""
    input:
        fasta="resources/phix/phix.fasta",
    output:
        expand("resources/phix/phix.{index}.bt2", index=range(1, 5)),
    log:
        "resources/phix/bowtie_build.log",
    params:
        prefix=lambda w, input: os.path.splitext(input.fasta)[0],
    threads: 1
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "bowtie2/2.4.5",
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.fasta} {params.prefix} >{log} 2>&1
        """


rule filter_phix_pe:
    """Maps reads against the phiX genome, keeping non-concordantly mapped"""
    input:
        bt_index=expand("resources/phix/phix.{index}.bt2", index=range(1, 5)),
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["phixfilt"]
        + ".fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["phixfilt"]
        + ".fastq.gz",
    output:
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["phixfilt"]
        + ".phixfilt.fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["phixfilt"]
        + ".phixfilt.fastq.gz",
    log:
        results
        + "/intermediate/preprocess/{sample}_{unit}_PHIX_pe"
        + preprocess_suffices["phixfilt"]
        + ".log",
    params:
        tmp_out=temppath,
        setting=config["bowtie2"]["extra_settings"],
        prefix="resources/phix/phix",
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=120,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "bowtie2/2.4.5",
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
            {params.tmp_out}/{wildcards.sample}_{wildcards.unit}_R%.filtered.fastq.gz > /dev/null 2>{log}
        mv {params.tmp_out}/{wildcards.sample}_{wildcards.unit}_R1.filtered.fastq.gz {output.R1}
        mv {params.tmp_out}/{wildcards.sample}_{wildcards.unit}_R2.filtered.fastq.gz {output.R2}
        """


rule filter_phix_se:
    """Maps reads against the phiX genome, keeping non-concordantly mapped"""
    input:
        bt_index=expand("resources/phix/phix.{index}.bt2", index=range(1, 5)),
        se=results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["phixfilt"]
        + ".fastq.gz",
    output:
        se=results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["phixfilt"]
        + ".phixfilt.fastq.gz",
    log:
        results
        + "/intermediate/preprocess/{sample}_{unit}_PHIX_se"
        + preprocess_suffices["phixfilt"]
        + ".log",
    params:
        tmp_out=temppath,
        setting=config["bowtie2"]["extra_settings"],
        prefix="resources/phix/phix",
    threads: config["bowtie2"]["threads"]
    resources:
        runtime=60,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "bowtie2/2.4.5",
    shell:
        """
        mkdir -p {params.tmp_out}
        bowtie2 \
            {params.setting} \
            -p {threads} \
            -x {params.prefix} \
            --un-gz \
            {params.tmp_out}/{wildcards.sample}_{wildcards.unit}_se.filtered.fastq.gz {input.se} > /dev/null 2>{log}
        mv {params.tmp_out}/{wildcards.sample}_{wildcards.unit}_se.filtered.fastq.gz {output.se}
        """


rule fastuniq:
    input:
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["fastuniq"]
        + ".fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["fastuniq"]
        + ".fastq.gz",
    output:
        R1=results
        + "/intermediate/preprocess/{sample}_{unit}_R1"
        + preprocess_suffices["fastuniq"]
        + ".fastuniq.fastq.gz",
        R2=results
        + "/intermediate/preprocess/{sample}_{unit}_R2"
        + preprocess_suffices["fastuniq"]
        + ".fastuniq.fastq.gz",
    log:
        results + "/intermediate/preprocess/{sample}_{unit}.fastuniq_pe.log",
    params:
        R1_intmp=temppath
        + "/{sample}_{unit}_R1"
        + preprocess_suffices["fastuniq"]
        + ".fastq",
        R2_intmp=temppath
        + "/{sample}_{unit}_R2"
        + preprocess_suffices["fastuniq"]
        + ".fastq",
        R1_outtmp=temppath
        + "/{sample}_{unit}_R1"
        + preprocess_suffices["fastuniq"]
        + ".fastuniq.fastq",
        R2_outtmp=temppath
        + "/{sample}_{unit}_R2"
        + preprocess_suffices["fastuniq"]
        + ".fastuniq.fastq",
        file_list=temppath + "/{sample}_{unit}.filelist",
    threads: 4
    resources:
        runtime=240,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    shell:
        """
        gunzip -c {input.R1} > {params.R1_intmp}
        gunzip -c {input.R2} > {params.R2_intmp}
        echo {params.R1_intmp} > {params.file_list}
        echo {params.R2_intmp} >> {params.file_list}
        fastuniq -i {params.file_list} -t q -o {params.R1_outtmp} \
            -p {params.R2_outtmp} >{log} 2>&1
        gzip -c {params.R1_outtmp} > {output.R1}
        gzip -c {params.R2_outtmp} > {output.R2}
        rm {params.R1_intmp} {params.R1_outtmp} {params.R2_intmp} {params.R2_outtmp}
        rm {params.file_list}
        """


rule fastuniq_se:
    """Dummy rule for fastuniq on single-end input"""
    input:
        se=results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["fastuniq"]
        + ".fastq.gz",
    output:
        se=results
        + "/intermediate/preprocess/{sample}_{unit}_se"
        + preprocess_suffices["fastuniq"]
        + ".fastuniq.fastq.gz",
    log:
        results + "/intermediate/preprocess/{sample}_{unit}.fastuniq_se.log",
    run:
        link(input.se, output.se)


rule fastqc:
    """Run fastqc on preprocessed data"""
    input:
        fastq=results
        + "/intermediate/preprocess/{sample}_{unit}_{pair}"
        + PREPROCESS
        + ".fastq.gz",
    output:
        zip=results
        + "/intermediate/fastqc/{sample}_{unit}_{pair}"
        + PREPROCESS
        + "_fastqc.zip",
    log:
        results + "/intermediate/fastqc/{sample}_{unit}_{pair}.log",
    params:
        dir=lambda w, output: os.path.dirname(output.zip),
    shadow:
        "shallow"
    resources:
        runtime=60,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
        if config["slurm_account"]
        else None,
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "FastQC/0.11.9",
    shell:
        """
        fastqc -q --noextract -o {params.dir} {input} >{log} 2>&1
        """


rule multiqc:
    """Summarize sample QC statistics in a report """
    input:
        multiqc_input(samples, config),
    output:
        html=report(
            results + "/report/samples_report.html",
            caption="../report/preprocessing.rst",
            category="Preprocessing",
        ),
        txt=results + "/report/samples_report_data/multiqc_general_stats.txt",
    log:
        results + "/report/multiqc.log",
    shadow:
        "shallow"
    params:
        config="workflow/.cfg/multiqc_preprocess_config.yaml",
        output_dir=lambda w, output: os.path.dirname(output.html),
    conda:
        "../envs/preprocess.yml"
    envmodules:
        "bioinfo-tools",
        "MultiQC/1.8",
    shell:
        """
        multiqc -f -c {params.config} -n samples_report.html \
            -o {params.output_dir} {input} >{log} 2>{log}
        """
