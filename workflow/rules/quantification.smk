localrules:
    write_featurefile,
    normalize_featurecount,
    aggregate_featurecount,
    sum_to_taxa,
    quantify_features,
    sum_to_rgi

rule write_featurefile:
    input:
        opj(config["results_path"],"annotation","{group}",
            "final_contigs.gff")
    output:
        opj(config["results_path"],"annotation","{group}",
            "final_contigs.features.gff")
    params:
        src=opj("source","utils","parse_prodigal_gff.py")
    shell:
        """
        python {params.src} {input} > {output}
        """

rule featurecount_pe:
    input:
        gff=opj(config["results_path"],"annotation","{group}",
                "final_contigs.features.gff"),
        bam=opj(config["results_path"],"assembly","{group}",
                "mapping","{sample}_{run}_pe"+POSTPROCESS+".bam")
    output:
        opj(config["results_path"],"assembly","{group}","mapping",
            "{sample}_{run}_pe.fc.tab"),
        opj(config["results_path"],"assembly","{group}","mapping",
            "{sample}_{run}_pe.fc.tab.summary")
    threads: 4
    params: tmpdir=config["tmpdir"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
    conda:
        "../../../envs/quantify.yml"
    shell:
        """
        featureCounts \
            -a {input.gff} \
            -o {output[0]} \
            -t CDS \
            -g gene_id \
            -M \
            -p \
            -B \
            -T {threads} \
            --tmpDir {params.tmpdir} \
            {input.bam}
        """

rule featurecount_se:
    input:
        gff=opj(config["results_path"],"annotation","{group}",
                "final_contigs.features.gff"),
        bam=opj(config["results_path"],"assembly","{group}",
                "mapping","{sample}_{run}_se"+POSTPROCESS+".bam")
    output:
        opj(config["results_path"],"assembly","{group}",
            "mapping","{sample}_{run}_se.fc.tab"),
        opj(config["results_path"],"assembly","{group}",
            "mapping","{sample}_{run}_se.fc.tab.summary")
    threads: 4
    params: tmpdir=config["tmpdir"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*30
    conda:
        "../../../envs/quantify.yml"
    shell:
        """
        featureCounts \
            -a {input.gff} \
            -o {output[0]} \
            -t CDS \
            -g gene_id \
            -M \
            -T {threads} \
            --tmpDir {params.tmpdir} \
            {input.bam}
        """

rule samtools_stats:
    input:
        opj(config["results_path"],"assembly","{group}",
                "mapping","{sample}_{run}_{seq_type}"+POSTPROCESS+".bam")
    output:
        opj(config["results_path"],"assembly","{group}",
                "mapping","{sample}_{run}_{seq_type}"+POSTPROCESS+".bam.stats")
    conda:
        "../../../envs/quantify.yml"
    shell:
        """
        samtools stats {input} > {output}
        """

rule normalize_featurecount:
    input:
        opj(config["results_path"],"assembly","{group}","mapping",
            "{sample}_{run}_{seq_type}.fc.tab"),
        opj(config["results_path"],"assembly","{group}","mapping",
            "{sample}_{run}_{seq_type}"+POSTPROCESS+".bam.stats")
    output:
        opj(config["results_path"],"assembly","{group}","mapping",
            "{sample}_{run}_{seq_type}.fc.tpm.tab"),
        opj(config["results_path"],"assembly","{group}","mapping",
            "{sample}_{run}_{seq_type}.fc.raw.tab")
    params:
        s="{sample}_{run}",
        src="source/utils/featureCountsTPM.py"
    shell:
        """
        rl=$(grep -w "average length" {input[1]} |  egrep -o "[0-9]+")
        python {params.src} \
            --rl $rl \
            -i {input[0]} \
            -o {output[0]} \
            --rc {output[1]} \
            --sampleName {params.s}
        """

rule aggregate_featurecount:
    """Aggregates feature count files and performs TPM normalization"""
    input:
        raw_files=get_all_files(samples, opj(config["results_path"], "assembly", "{group}", "mapping"), ".fc.raw.tab"),
        tpm_files=get_all_files(samples, opj(config["results_path"], "assembly", "{group}", "mapping"), ".fc.tpm.tab"),
        gff_file=opj(config["results_path"],"annotation","{group}","final_contigs.features.gff")
    output:
        raw=opj(config["results_path"],"annotation","{group}","fc.raw.tab"),
        tpm=opj(config["results_path"],"annotation","{group}","fc.tpm.tab")
    run:
        gff_df=pd.read_csv(input.gff_file, header=None, usecols=[0,8], names=["contig","gene"], sep="\t")
        gff_df=gff_df.assign(gene_id=pd.Series([x.replace("gene_id ","") for x in gff_df.gene], index=gff_df.index))
        gff_df=gff_df.assign(suffix=pd.Series([x.split(" ")[-1].split("_")[-1] for x in gff_df.gene],index=gff_df.index))
        gff_df=gff_df.assign(orf=pd.Series(gff_df.contig+"_"+gff_df.suffix, index=gff_df.index))
        gff_df=gff_df[["orf","gene_id"]]

        raw_df=concat_files(input.raw_files, gff_df)
        tpm_df=concat_files(input.tpm_files, gff_df)
        raw_df.to_csv(output.raw, sep="\t")
        tpm_df.to_csv(output.tpm, sep="\t")

rule quantify_features:
    input:
        abund=opj(config["results_path"],"annotation","{group}","fc.{fc_type}.tab"),
        annot=opj(config["results_path"],"annotation","{group}","{db}.parsed.tab")
    output:
        opj(config["results_path"],"annotation","{group}","{db}.parsed.{fc_type}.tab")
    shell:
        """
        python source/utils/eggnog-parser.py \
            quantify {input.abund} {input.annot} {output[0]}
        """

rule sum_to_taxa:
    input:
        tax=opj(config["results_path"],"annotation","{group}","taxonomy",
            "orfs.{db}.taxonomy.tsv".format(db=config["taxdb"])),
        count=opj(config["results_path"],"annotation","{group}","fc.raw.tab"),
        norm=opj(config["results_path"],"annotation","{group}","fc.tpm.tab")
    output:
        count=opj(config["results_path"],"annotation","{group}","taxonomy","tax.raw.tab"),
        norm=opj(config["results_path"],"annotation","{group}","taxonomy","tax.tpm.tab")
    run:
        header=["protein","superkingdom", "phylum","class","order","family","genus","species"]
        df=pd.read_csv(input.tax, sep="\t", index_col=0, header=None, names=header)
        count_df=pd.read_csv(input.count, header=0, index_col=0, sep="\t")
        norm_df=pd.read_csv(input.norm, header=0, index_col=0, sep="\t")

        taxa_count=pd.merge(df,count_df,right_index=True,left_index=True)
        taxa_count_sum=taxa_count.groupby(header[1:]).sum().reset_index()
        taxa_count_sum.to_csv(output.count, sep="\t", index=False)

        taxa_norm=pd.merge(df,norm_df,right_index=True,left_index=True)
        taxa_norm_sum=taxa_norm.groupby(header[1:]).sum().reset_index()
        taxa_norm_sum.to_csv(output.norm, sep="\t", index=False)

rule sum_to_rgi:
    input:
        annot=opj(config["results_path"], "annotation", "{group}", "rgi.out.txt"),
        abund=opj(config["results_path"],"annotation","{group}","fc.{fc_type}.tab")
    output:
        opj(config["results_path"], "annotation", "{group}", "rgi.{fc_type}.tab")
    run:
        annot=pd.read_csv(input.annot, sep="\t", header=0, index_col=0, usecols=[0,16])
        # Rename index for annotations to remove text after whitespace
        annot.rename(index=lambda x: x.split(" ")[0], inplace=True)
        abund=pd.read_csv(input.abund, sep="\t", header=0, index_col=0)
        df=pd.merge(annot, abund, left_index=True, right_index=True)
        # Sum to Gene family
        dfsum=df.groupby("AMR Gene Family").sum()
        dfsum.to_csv(output[0], sep="\t", index=True, header=True)
