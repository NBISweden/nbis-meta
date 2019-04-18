binning_input = []

if config["maxbin"]:
    binning_input += expand(opj(config["results_path"],"maxbin","{group}","{l}","summary_stats.tsv"),
        group = assemblyGroups.keys(), l=config["min_contig_length"])
if config["concoct"]:
    binning_input += expand(opj(config["results_path"],"concoct","{group}","{l}","summary_stats.tsv"),
        group = assemblyGroups.keys(), l=config["min_contig_length"])
if config["metabat"]:
    binning_input += expand(opj(config["results_path"],"metabat","{group}","{l}","summary_stats.tsv"),
        group = assemblyGroups.keys(),
        l=[x for x in config["min_contig_length"] if x>=1500]) # Metabat has a minlength of 1500

include: "../rules/Binning/coverage.rules"
include: "../rules/Binning/binning.rules"