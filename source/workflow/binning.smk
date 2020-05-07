binning_input = []

if config["maxbin"]:
    binning_input += expand(opj(config["results_path"],
                                "binning", "maxbin","{group}","{l}","summary_stats.tsv"),
                            group=assemblyGroups.keys(),
                            l=config["min_contig_length"])
if config["concoct"] and platform.uname().system != "Darwin":
    binning_input += expand(opj(config["results_path"],
                                "binning", "concoct","{group}","{l}","summary_stats.tsv"),
                            group=assemblyGroups.keys(),
                            l=config["min_contig_length"])
if config["metabat"]:
    binning_input += expand(opj(config["results_path"],
                                "binning", "metabat","{group}","{l}","summary_stats.tsv"),
                            group=assemblyGroups.keys(),
                            l=config["min_contig_length"]) # Metabat has a minlength of 1500

if config["checkm"]:
    binning_input.append(opj(config["report_path"], "checkm", "checkm.stats.tsv"))
    binning_input.append(opj(config["report_path"], "checkm", "checkm.profiles.tsv"))
    include: "../rules/Binning/checkm.smk"

if config["gtdbtk"] and platform.uname().system != "Darwin":
    binning_input.append(opj(config["report_path"], "gtdbtk", "gtdbtk.summary.tsv"))
    binning_input.append(opj(config["report_path"], "mag_annotation", "tRNA.total.tsv"))
    binning_input.append(opj(config["report_path"], "mag_annotation", "rRNA.types.tsv"))
    include: "../rules/Binning/gtdbtk.smk"
    include: "../rules/Binning/annotate_bins.smk"

include: "../rules/Binning/coverage.smk"
include: "../rules/Binning/binning.smk"