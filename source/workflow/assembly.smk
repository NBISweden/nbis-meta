assembly_input = expand(opj(config["results_path"],"assembly","{group}",
                   "final_contigs.fa"),group=assemblyGroups.keys())
assembly_input += expand(opj(config["results_path"],"assembly","{group}",
                   "final_contigs.bed"),group=assemblyGroups.keys())
assembly_input += [opj(config["report_path"],"assembly_stats.txt"),
                 opj(config["report_path"],"assembly_size_dist.txt")]

include: "../rules/Assembly/assembly.smk"

if not config["annotation"]:
    include: "../rules/Annotation/orfcalling.rules"
    include: "../rules/Annotation/markduplicates.smk"
    include: "../rules/Annotation/quantification.rules"