# Workflow module for protein annotations
annotation_input = []

for group in assemblyGroups.keys():
    # Add orfcalling results
    annotation_input.append(opj(config["results_path"],"annotation",
                                group,"final_contigs.gff"))
    if config["infernal"]:
        annotation_input.append(opj(config["results_path"],"annotation",
                                    group,"final_contigs.cmscan"))
    if config["tRNAscan"]:
        annotation_input.append(opj(config["results_path"],"annotation",
                                    group,"tRNA.out"))
    # Add EGGNOG annotation
    if config["eggnog"]:
        annotation_input+=expand(opj(config["results_path"],"annotation",
                                     group,"{db}.parsed.{fc}.tab"),
                                 db=["enzymes","pathways","kos","modules"],
                                 fc=["raw","tpm"])
    # Add PFAM annotation
    if config["pfam"]:
        annotation_input+=expand(opj(config["results_path"],"annotation",
                                    group,"pfam.parsed.{fc}.tab"),
                                fc=["tpm","raw"])
    # Add taxonomic annotation
    if config["taxonomic_annotation"]:
        annotation_input+=expand(opj(config["results_path"],"annotation",
                                     group,"taxonomy","tax.{fc}.tab"),
                                 fc=["tpm","raw"])
    # Add Resistance Gene Identifier output
    if config["rgi"]:
        annotation_input += expand(opj(config["results_path"],"annotation",
                                       group,"rgi.{fc}.tab"),
                                   fc=["raw","tpm"])
        annotation_input.append(opj(config["results_path"],"annotation",
                                    group,"rgi.out.txt"))

include: "../rules/Annotation/markduplicates.smk"
include: "../rules/Annotation/orfcalling.smk"
include: "../rules/Annotation/prot_annotation.smk"
include: "../rules/Annotation/quantification.smk"
include: "../rules/Annotation/taxonomy.smk"