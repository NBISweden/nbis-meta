localrules:
    classifier2krona,
    all2krona

rule classifier2krona:
    input:
        opj(config["results_path"],"{classifier}",
            "{sample}_{run}_{seq_type}.kreport"),
        opj("resources","krona","taxonomy.tab")
    output:
        opj(config["results_path"],"{classifier}",
            "{sample}_{run}_{seq_type}.html")
    params:
        tax="resources/krona"
    conda:
        "../../../envs/krona.yml"
    shell:
        """
        ktImportTaxonomy -t 5 -m 3 \
            -tax {params.tax} -o {output[0]} \
            {input[0]},{wildcards.sample}_{wildcards.run}
        """

def get_krona_input(samples, classifier):
    input_string=""
    files = get_all_files(samples,opj(config["results_path"],
                                    classifier),".kreport")
    for f in files:
        sample_run=os.path.basename(f).replace("_pe.kreport","").replace("_se.kreport","")
        input_string+=" {},{}".format(f,sample_run)
    return input_string


rule all2krona:
    input:
        f=get_all_files(samples,opj(config["results_path"],
                                    "{classifier}"),".kreport"),
        h=get_all_files(samples,opj(config["results_path"],
                                    "{classifier}"),".html"),
        t=opj("resources","krona","taxonomy.tab")
    output:
        opj(config["report_path"],"{classifier}","{classifier}.krona.html")
    params:
        tax="resources/krona",
        input_string=get_krona_input(samples, "{classifier}")
    conda:
        "../../../envs/krona.yml"
    shell:
         """
         ktImportTaxonomy \
            -t 5 -m 3 -tax {params.tax} -o {output[0]} {params.input_string}
         """