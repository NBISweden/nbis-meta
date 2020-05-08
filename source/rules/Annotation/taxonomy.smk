localrules: tango_assign_orfs, download_sourmash_db, sourmash_compute, merge_tango_sourmash

rule tango_search:
    input:
        db=opj(config["resource_path"],config["taxdb"],"diamond.dmnd"),
        fasta=opj(config["results_path"],"assembly","{group}",
                  "final_contigs.fa")
    output:
        opj(config["results_path"],"annotation","{group}",
            "final_contigs.{db}.tsv.gz".format(db=config["taxdb"]))
    log:
        opj(config["results_path"],"annotation","{group}", "tango_search.log")
    params:
        tmpdir = config["scratch_path"],
        min_len = config["taxonomy_min_len"],
        settings = config["tango_search_params"]
    threads: config["diamond_threads"]
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*10
    conda:
        "../../../envs/tango.yml"
    shell:
        """
        tango search {params.settings} -p {threads} \
            --tmpdir {params.tmpdir} -l {params.min_len} \
            {input.fasta} {input.db} {output[0]} >{log} 2>&1
        """

rule tango_assign:
    input:
        opj(config["results_path"],"annotation","{group}",
            "final_contigs.{db}.tsv.gz".format(db=config["taxdb"])),
        opj(config["resource_path"],"taxonomy","taxonomy.sqlite")
    output:
        opj(config["results_path"],"annotation","{group}","taxonomy",
            "tango.{db}.taxonomy.tsv".format(db=config["taxdb"]))
    log:
        opj(config["results_path"],"annotation","{group}","taxonomy",
            "tango_assign.log")
    params:
        taxonomy_ranks = config["taxonomy_ranks"],
        taxdir = opj(config["resource_path"],"taxonomy"),
        settings = config["tango_assign_params"]
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*6
    conda:
        "../../../envs/tango.yml"
    shell:
         """
         tango assign {params.settings} -p {threads} -m rank_lca \
            --reportranks {params.taxonomy_ranks} -t {params.taxdir} \
            {input[0]} {output[0]} > {log} 2>&1
         """

## Sourmash rules

rule sourmash_compute:
    input:
        opj(config["results_path"],"assembly","{group}", "final_contigs.fa")
    output:
        opj(config["results_path"],"assembly","{group}", "final_contigs.fa.sig")
    log:
        opj(config["results_path"],"assembly","{group}", "sourmash_compute.log")
    conda:
        "../../../envs/sourmash.yml"
    params:
        frac = config["sourmash_fraction"]
    shell:
        """
        sourmash compute --singleton --scaled {params.frac} \
            -k 31 -o {output} {input} > {log} 2>&1
        """

rule sourmash_classify:
    input:
        sig = opj(config["results_path"],"assembly","{group}",
                 "final_contigs.fa.sig"),
        db = opj(config["resource_path"], "sourmash", "genbank-k31.lca.json")
    output:
        csv = opj(config["results_path"], "annotation", "{group}", "taxonomy",
                  "sourmash.taxonomy.csv")
    log:
        opj(config["results_path"], "annotation", "{group}", "taxonomy",
            "sourmash.log")
    params:
        frac = config["sourmash_fraction"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*30
    conda:
        "../../../envs/sourmash.yml"
    shell:
        """
        sourmash lca classify --db {input.db} --scaled {params.frac} \
            --query {input.sig} -o {output.csv} > {log} 2>&1
        """

rule merge_tango_sourmash:
    input:
        smash = opj(config["results_path"], "annotation", "{group}",
                    "taxonomy", "sourmash.taxonomy.csv"),
        tango = opj(config["results_path"],"annotation","{group}",
                    "taxonomy", "tango.{db}.taxonomy.tsv".format(db=config["taxdb"]))
    output:
        tax = opj(config["results_path"],"annotation","{group}","taxonomy",
        "final_contigs.taxonomy.tsv")
    log:
        opj(config["results_path"],"annotation","{group}","taxonomy","merge.log")
    shell:
        """
        python source/utils/tango_mash.py \
            {input.smash} {input.tango} > {output.tax} 2>{log}
        """

rule tango_assign_orfs:
    input:
        tax=opj(config["results_path"],"annotation","{group}","taxonomy",
            "final_contigs.taxonomy.tsv"),
        gff=opj(config["results_path"],"annotation","{group}",
                "final_contigs.gff")
    output:
        tax=opj(config["results_path"],"annotation","{group}","taxonomy",
            "orfs.{db}.taxonomy.tsv".format(db=config["taxdb"]))
    run:
        import pandas as pd
        gff_df=pd.read_csv(input.gff, header=None, sep="\t", comment="#",
                           usecols=[0,8], names=["contig","id"])
        # Extract ids
        ids=["{}_{}".format(gff_df.loc[i, "contig"], gff_df.loc[i, "id"].split(";")[0].split("_")[-1]) for i in gff_df.index]
        gff_df.loc[:, "id"]=ids
        # Read taxonomy for contigs
        tax_df=pd.read_csv(input.tax, header=0, sep="\t", index_col=0)
        # Merge dataframes
        orf_tax_df=pd.merge(gff_df, tax_df, left_on="contig",
                            right_index=True, how="outer")
        # When using 'outer' merging there may be contigs with no called ORF
        # but with a tax assignment. Drop these contigs.
        orf_tax_df = orf_tax_df.loc[orf_tax_df["id"]==orf_tax_df["id"]]
        # Set Unclassified for NA values
        orf_tax_df.fillna("Unclassified", inplace=True)
        # Set index to ORF ids
        orf_tax_df.set_index("id", inplace=True)
        orf_tax_df.drop("contig", axis=1, inplace=True)
        orf_tax_df.to_csv(output.tax, sep="\t", index=True, header=True)