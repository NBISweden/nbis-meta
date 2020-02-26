localrules:
    write_bed,
    plot_assembly_stats,
    stat_alignment_rate,
    plot_group_mean_alignment_rate,
    assembly_stats

#############
# Functions #
#############

def get_all_group_files(g):
  files = []
  for sample in assemblyGroups[g].keys():
    for run in assemblyGroups[g][sample].keys():
      for pair in assemblyGroups[g][sample][run].keys():
        files.append(assemblyGroups[g][sample][run][pair][0])
  return files

def get_paired_group_files(g, r):
    files = []
    for sample in assemblyGroups[g].keys():
        for run in assemblyGroups[g][sample].keys():
            if not is_pe(assemblyGroups[g][sample][run]):
                continue
            files.append(assemblyGroups[g][sample][run][r])
    return sorted(files)

def get_bamfiles(g):
  files = []
  for sample in assemblyGroups[g].keys():
    for run in assemblyGroups[g][sample].keys():
      if "R2" in assemblyGroups[g][sample][run].keys():
        files.append(opj(config["results_path"],"assembly",g,"mapping",sample+"_"+run+"_pe"+POSTPROCESS+".bam"))
      else:
        files.append(opj(config["results_path"],"assembly",g,"mapping",sample+"_"+run+"_se"+POSTPROCESS+".bam"))
  return files

def add_series(df,html):
    html+= "\tseries: [\n"
    for assembly in df.assembly.unique():
        r = df.loc[df.assembly==assembly]
        html+="\t\t{{\n\t\tname: '{}',\n".format(assembly)
        html+="\t\tdata: ["
        for i in r.index:
            min_length,percent = r.loc[i,["min_length","%"]]
            html+="[{},{}],".format(min_length,percent)
        html = html.rstrip(",")
        html+="]\n\t\t},\n"
    html+="\t\t],\n"
    return html

def generate_html(df):
  with open("source/reporting/highcharts_template.html") as fh:
    template = fh.read()
    html = ""
    for line in template.split("\n"):
      if line=="#ADD_SERIES":
        html = add_series(df,html)
      else:
        html+="{}\n".format(line)
  return html

############
# Assembly #
############

def rename_records(f, fh, i):
    """
    Prepends a number to read ids

    :param f: Input fastq file (gzipped)
    :param fh: Output filehandle
    :param i: File index to prepend to reads
    :return: Output filehandle
    """
    for record in SeqIO.parse(gz.open(f, 'rt'), 'fastq'):
        record.id = "{}_{}".format(i, record.id)
        SeqIO.write(record, fh, "fastq")
    return fh

if config["metaspades"]:
    rule generate_metaspades_input_list:
        """Generate input files for use with Metaspades"""
        input:
            lambda wildcards: get_all_group_files(wildcards.group)
        output:
            R1=temp(opj(config["results_path"],"assembly",
                        "{group}","R1.fq")),
            R2=temp(opj(config["results_path"],"assembly",
                        "{group}","R2.fq")),
            se=touch(temp(opj(config["results_path"],"assembly",
                        "{group}","se.fq")))
        run:
            files = {"R1": [], "R2": [], "se": []}
            # Collect all files belonging to the assembly group
            for sample in assemblyGroups[wildcards.group].keys():
                for run in assemblyGroups[wildcards.group][sample]:
                    for pair in assemblyGroups[wildcards.group][sample][run].keys():
                        files[pair].append(assemblyGroups[wildcards.group][sample][run][pair][0])
            # Rename and concatenate reads (required for Metaspades)
            with open(output.R1, 'w') as fh1, open(output.R2, 'w') as fh2, open(output.se, 'w') as fhse:
                for i, f in enumerate(files["R1"]):
                    f2=files["R2"][i]
                    fh1 = rename_records(f, fh1, i)
                    fh2 = rename_records(f2, fh2, i)
                for i, f in enumerate(files["se"], start=i+1):
                    fhse = rename_records(f, fhse, i)

    rule run_metaspades:
        input:
            R1=opj(config["results_path"],"assembly",
                        "{group}","R1.fq"),
            R2=opj(config["results_path"],"assembly",
                        "{group}","R2.fq"),
            se=opj(config["results_path"],"assembly",
                        "{group}","se.fq")
        output:
            opj(config["results_path"],"assembly","{group}","final_contigs.fa"),
            opj(config["results_path"],"assembly","{group}","spades.log")
        params:
            intermediate_contigs=opj(config["intermediate_path"],"assembly",
                                     "{group}","intermediate_contigs"),
            corrected=opj(config["intermediate_path"],"assembly",
                          "{group}","corrected"),
            additional_settings=config["metaspades_additional_settings"],
            tmp=opj(config["tmpdir"],"{group}.metaspades"),
            output_dir=opj(config["results_path"],"assembly","{group}")
        threads: config["assembly_threads"]
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60*4
        conda:
            "../../../envs/metaspades.yml"
        shell:
            """
            # Create directories
            mkdir -p {params.tmp}
            # Only use single-end if present
            if [ -s {input.se} ]; then
                single="-s {input.se}"
            else
                single=""
            fi
            metaspades.py \
                -t {threads} \
                -1 {input.R1} \
                -2 {input.R2} \
                $single \
                -o {params.tmp}
            
            # If set to keep intermediate contigs, move to intermediate folder before deleting
            if [ "{config[metaspades_keep_intermediate]}" == "True" ]; then
                mkdir -p {params.intermediate_contigs}
                rsync -azv {params.tmp}/K* {params.intermediate_contigs}
            fi
            if [ "{config[metaspades_keep_corrected]}" == "True" ]; then
                mkdir -p {params.corrected}
                rsync -azv {params.tmp}/corrected {params.corrected}
            fi
            
            # Clear intermediate contigs
            rm -rf {params.tmp}/K*
            # Clear corrected reads dir
            rm -rf {params.tmp}/corrected
            # Sync tmp output to outdir before removing
            rsync -azv {params.tmp}/* {params.output_dir}
            rm -rf {params.tmp}
            mv {params.output_dir}/scaffolds.fasta {params.output_dir}/final_contigs.fa
            """


else:
    rule generate_megahit_input_list:
        """Generate input lists for Megahit"""
        input:
            lambda wildcards: get_all_group_files(wildcards.group)
        output:
            R1=temp(opj(config["results_path"],"assembly",
                            "{group}","input_1")),
            R2=temp(opj(config["results_path"],"assembly",
                            "{group}","input_2")),
            se=temp(opj(config["results_path"],"assembly",
                            "{group}","input_se"))
        run:
            files = {"R1": [], "R2": [], "se": []}
            for sample in assemblyGroups[wildcards.group].keys():
                for run in assemblyGroups[wildcards.group][sample]:
                    for pair in assemblyGroups[wildcards.group][sample][run].keys():
                        files[pair].append(assemblyGroups[wildcards.group][sample][run][pair][0])
            with open(output.R1, 'w') as fh1, open(output.R2, 'w') as fh2, open(output.se, 'w') as fhse:
                fh1.write(",".join(files["R1"]))
                fh2.write(",".join(files["R2"]))
                fhse.write(",".join(files["se"]))

    rule run_megahit:
        input:
            R1=opj(config["results_path"],"assembly",
                            "{group}","input_1"),
            R2=opj(config["results_path"],"assembly",
                            "{group}","input_2"),
            se=opj(config["results_path"],"assembly",
                            "{group}","input_se")
        output:
            opj(config["results_path"],"assembly","{group}","final_contigs.fa"),
            opj(config["results_path"],"assembly","{group}","log")
        params:
            intermediate_contigs=opj(config["intermediate_path"],"assembly",
                                     "{group}","intermediate_contigs"),
            additional_settings=config["megahit_additional_settings"],
            tmp=opj(config["tmpdir"],"{group}.megahit"),
            output_dir=opj(config["results_path"],"assembly","{group}")
        threads: config["assembly_threads"]
        resources:
            runtime=lambda wildcards, attempt: attempt**2*60*4
        conda:
            "../../../envs/megahit.yml"
        shell:
            """
            rm -rf {params.tmp}
            
            # Get input strings
            R1=$(cat {input.R1})
            R2=$(cat {input.R2})
            # Only use single-end if present
            se=$(cat {input.se})
            if [ -s {input.se} ]; then
                single="-r $se"
            else
                single=""
            fi
            
            # Run Megahit
            megahit \
                -t {threads} \
                -1 $R1 \
                -2 $R2 \
                $single \
                -o {params.tmp} \
                {params.additional_settings} \
                >/dev/null 2>&1
            
            # Sync intermediate contigs if asked for
            if [ "{config[megahit_keep_intermediate]}" == "True" ]; then
                mkdir -p {params.intermediate_contigs}
                rsync -azv {params.tmp}/intermediate_contigs/* {params.intermediate_contigs}
            fi
            
            # Cleanup intermediate
            rm -rf {params.tmp}/intermediate_contigs
            
            # Sync tmp output to outdir before removing
            rsync -azv {params.tmp}/* {params.output_dir}
            rm -rf {params.tmp}
            mv {params.output_dir}/final.contigs.fa {params.output_dir}/final_contigs.fa
            """

rule write_bed:
    """Creates bed-format file from assembly"""
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        opj(config["results_path"],"assembly","{group}","final_contigs.bed")
    shell:
        """
        python source/utils/fasta2bed.py -i {input} -o {output}
        """

###########
# Mapping #
###########

rule bowtie_build:
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        expand(opj(config["results_path"],"assembly","{{group}}","final_contigs.fa.{index}.bt2l"),index=range(1,5))
    params: prefix = opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    threads: config["bowtie2_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        bowtie2-build --large-index --threads {threads} {params.prefix} {params.prefix}
        """

rule bowtie_map_pe:
    input:
        bt_index = expand(opj(config["results_path"],"assembly","{{group}}","final_contigs.fa.{index}.bt2l"),index=range(1,5)),
        R1 = expand(opj(config["intermediate_path"], "preprocess", "{{sample}}_{{run}}_R1{p}.fastq.gz"), p=PREPROCESS),
        R2 = expand(opj(config["intermediate_path"], "preprocess", "{{sample}}_{{run}}_R2{p}.fastq.gz"), p=PREPROCESS)
    output:
        bam=temp(opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_pe.bam")),
        log=opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_pe.bam.log")
    params:
        temp_bam=opj(config["tmpdir"],"{group}-mapping-{sample}_{run}_pe.bam"),
        setting = config["bowtie2_params"],
        prefix = opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    threads: config["bowtie2_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -1 {input.R1} -2 {input.R2} 2> {output.log} | \
        samtools view -bh - | samtools sort -o {params.temp_bam}
        mv {params.temp_bam} {output.bam}
        """

rule bowtie_map_se:
    input:
        bt_index = expand(opj(config["results_path"],"assembly","{{group}}","final_contigs.fa.{index}.bt2l"),index=range(1,5)),
        se = expand(opj(config["intermediate_path"], "preprocess", "{{sample}}_{{run}}_se{p}.fastq.gz"), p=PREPROCESS)
    output:
        bam=temp(opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_se.bam")),
        log=opj(config["results_path"],"assembly","{group}","mapping","{sample}_{run}_se.bam.log")
    params:
        temp_bam=opj(config["tmpdir"],"{group}-mapping-{sample}_{run}_se.bam"),
        setting = config["bowtie2_params"],
        prefix = opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    threads: config["bowtie2_threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    shell:
        """
        bowtie2 {params.setting} -p {threads} -x {params.prefix} -U {input.se} 2>{output.log} | samtools view -bh - \
         | samtools sort -o {params.temp_bam}
        mv {params.temp_bam} {output.bam}
        """

##############
# Statistics #
##############

rule assembly_stats:
    input:
        expand(opj(config["results_path"],"assembly","{group}","final_contigs.fa"), group = assemblyGroups.keys())
    output:
        opj(config["report_path"],"assembly_stats.txt"),
        opj(config["report_path"],"assembly_size_dist.txt")
    run:
        from source.utils import assembly_stats
        stat_result = pd.DataFrame()
        sizedist_result = pd.DataFrame()
        for f in input:
            print("Generating stats for {}".format(f))
            name = f.split("/")[-2]
            contig_lengths = assembly_stats.store_lengths(f)
            stat_df = assembly_stats.generate_stat_df(contig_lengths)
            size_dist = assembly_stats.size_distribute(contig_lengths)
            stat_df["assembly"] = [name]*len(stat_df)
            size_dist["assembly"] = [name]*len(size_dist)
            stat_result = pd.concat([stat_result,stat_df])
            sizedist_result = pd.concat([sizedist_result,size_dist])
        stat_result = stat_result[["assembly","contigs", "total_size_bp", "min_length", "max_length", "avg_length", "median_length",
             "N50_length", "N90_length"]]
        stat_result.to_csv(output[0],sep="\t",index=False)
        sizedist_result = sizedist_result[["assembly","min_length","num_contigs","total_length","%"]]
        sizedist_result.to_csv(output[1],sep="\t",index=False)

rule plot_assembly_stats:
    input:
        opj(config["report_path"],"assembly_stats.txt"),
        opj(config["report_path"],"assembly_size_dist.txt")
    output:
        opj(config["report_path"],"assembly_stats.pdf"),
        opj(config["report_path"],"assembly_size_dist.pdf"),
        opj(config["report_path"],"assembly_size_distribution.html")
    run:
        import matplotlib as mpl
        mpl.use('agg')
        import pandas as pd, seaborn as sns, matplotlib.pyplot as plt
        plt.style.use('ggplot')
        stat_result = pd.read_csv(input[0], sep="\t", header=0)
        sizedist_result = pd.read_csv(input[1], sep="\t", header=0)

        stat_result_m = pd.melt(stat_result, id_vars=["assembly"])
        ax = sns.factorplot(kind="bar", x="assembly", y="value", col="variable", data=stat_result_m, size=3,
                    sharey=False, col_wrap=4)
        ax.set_xticklabels(rotation=90)
        ax.set_titles("{col_name}")
        plt.savefig(output[0], dpi=300, bbox_inches="tight")

        ax = sns.factorplot(data=sizedist_result, hue="assembly", x="min_length", y="%", size=4, scale=0.6,aspect=1.2)
        ax.set_xticklabels(rotation=270)
        ax.set_ylabels("% of assembly")
        ax.axes[0][0].set_ylim(-5,105);
        plt.savefig(output[1],dpi=300, bbox_inches="tight")

        html = generate_html(sizedist_result)
        with open(output[2], 'w') as fh:
          fh.write(html)

rule stat_alignment_rate:
    input:
        lambda wildcards: get_bamfiles(wildcards.group)
    output:
        temp(opj(config["report_path"],"{group}_alignment_rate.tab"))
    run:
        import re
        ali_re = re.compile(" mapped \((\d+.\d+)%")
        result = {}
        for f in input:
            for line in shell("samtools flagstat {f}", iterable = True):
                try:
                    ali = ali_re.search(line).groups()[0]
                    break
                except AttributeError:
                    continue
            basename = os.path.basename(f)
            sample_run = basename.replace("_pe{post}.bam".format(post=POSTPROCESS),"")
            sample_run = sample_run.replace("_se{post}.bam".format(post=POSTPROCESS),"")
            result[sample_run] = ali
        df = pd.DataFrame(result,index = ["alignment_rate"]).T
        df.index.name = wildcards.group
        df.to_csv(output[0],sep="\t")

rule plot_group_mean_alignment_rate:
    input:
        expand(opj(config["report_path"],"{group}_alignment_rate.tab"), group=assemblyGroups.keys())
    output:
        plot=opj(config["report_path"],"alignment_rate.pdf"),
        table=opj(config["report_path"],"alignment_rate.tab")
    run:
        import matplotlib as mpl
        mpl.use('agg')
        import pandas as pd, seaborn as sns, matplotlib.pyplot as plt
        plt.style.use('ggplot')
        df = pd.DataFrame()
        for f in input:
            _df = pd.read_csv(f,header=0,sep="\t")
            _df["assembly_group"] = [_df.columns[0]]*len(_df)
            _df.columns = ["sample_run","alignment_rate","assembly_group"]
            df = pd.concat([df,_df])
        value_name = "% overall alignment"
        dfm = pd.melt(df,id_vars=["sample_run","assembly_group"],value_vars=["alignment_rate"], value_name=value_name)
        fig,axes = plt.subplots(ncols=1, nrows=2, sharex=False, sharey=False, figsize=(6,12))
        ax1 = sns.barplot(data=dfm, x="assembly_group", y=value_name, errwidth=1, ax=axes[0])
        ax2 = sns.barplot(data=dfm, x="sample_run", y=value_name, hue="assembly_group", ax=axes[1])
        ax1.set_xticklabels(ax1.get_xticklabels(),rotation=270)
        ax1.set_title("Groupwise alignment rate");
        ax2.set_xticklabels(ax2.get_xticklabels(),rotation=270,fontdict={'size': 8})
        ax2.set_title("Samplewise alignment rate");
        ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(output.plot, dpi=300, bbox_inches="tight")
        dfm.to_csv(output.table, sep="\t")

rule assembly_all:
    input:
        expand(opj(config["results_path"],"assembly","{group}","final_contigs.fa"),group=assemblyGroups.keys()),
        expand(opj(config["results_path"],"assembly","{group}","final_contigs.bed"),group=assemblyGroups.keys()),
        opj(config["report_path"],"assembly_stats.pdf"),
        opj(config["report_path"],"assembly_size_dist.pdf"),
        opj(config["report_path"],"assembly_size_distribution.html"),
        opj(config["report_path"],"assembly_stats.txt"),
        opj(config["report_path"],"assembly_size_dist.txt"),
        opj(config["report_path"],"assembly_standard_err_cov.pdf")
    output:
        expand(opj(config["results_path"],"assembly","{group}","assembly.done"),group=assemblyGroups.keys())
    run:
        for f in output:
            shell("touch {f}")
