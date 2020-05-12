# Workflow for assembly
assembly_input = expand(opj(config["results_path"], "assembly", "{group}", "final_contigs.{suff}"),
                        group=assemblyGroups.keys(), suff = ["fa", "bed"])
assembly_input.append(opj(config["report_path"], "assembly", "assembly_stats.pdf"))
assembly_input.append(opj(config["report_path"], "assembly", "assembly_size_dist.pdf"))
assembly_input.append(opj(config["report_path"], "assembly", "alignment_frequency.pdf"))

include: "../rules/Assembly/assembly.smk"

if not config["annotation"]:
    include: "../rules/Annotation/orfcalling.smk"
    include: "../rules/Annotation/markduplicates.smk"
    include: "../rules/Annotation/quantification.smk"