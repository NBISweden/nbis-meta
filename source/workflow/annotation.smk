annotation_input = []

# if any assembly groups have been specified in the sample file then run these
for group in assemblyGroups.keys():
    # Add orfcalling results
    annotation_input.append(opj(config["results_path"],"annotation",group,"final_contigs.gff"))
    if config["infernal"]:
        annotation_input.append(opj(config["results_path"],"annotation",group,"final_contigs.cmscan"))
    if config["tRNAscan"]:
        annotation_input.append(opj(config["results_path"],"annotation",group,"tRNA.out"))
    # Add EGGNOG annotation
    if config["eggnog"]:
        for parsed in ["enzymes","pathways","kos","modules"]:
            annotation_input.append(opj(config["results_path"],"annotation",group,"{}.parsed.tpm.tab".format(parsed)))
        # Add normalized output for modules and pathways
        annotation_input += expand(opj(config["results_path"], "annotation", group,"{eggnog_type}.parsed.{fc_type}.normalized.tab"),
                                   eggnog_type = ["modules","pathways"], fc_type = ["count", "tpm"])
    # Add PFAM annotation
    if config["pfam"]:
        annotation_input.append(opj(config["results_path"],"annotation",group,"pfam.parsed.tpm.tab"))
    # Add taxonomic annotation
    if config["taxonomic_annotation"]:
        annotation_input.append(opj(config["results_path"],"annotation",group,"taxonomy","taxonomy.tpm.krona.html"))
        annotation_input.append(opj(config["results_path"],"annotation",group,"taxonomy","tax.tpm.tab"))
    # Add Resistance Gene Identifier output
    if config["rgi"]:
        annotation_input += expand(opj(config["results_path"],"annotation",group,"rgi.{fc_type}.tab"), fc_type = ["count","tpm"])
        annotation_input.append(opj(config["results_path"],"annotation",group,"rgi.out.txt"))

include: "../rules/Annotation/markduplicates.smk"
include: "../rules/Annotation/orfcalling.smk"
include: "../rules/Annotation/prot_annotation.rules"
include: "../rules/Annotation/quantification.rules"
include: "../rules/Annotation/taxonomic_annotation.rules"