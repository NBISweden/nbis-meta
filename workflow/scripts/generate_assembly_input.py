#!/usr/bin/env python


def metaspades_input(sm):
    """
    Generates fastq files to use as input for metaspades assembler

    :param sm: snakemake object
    :return:
    """
    from common import rename_records
    files = {"R1": [], "R2": [], "se": []}
    assembly_dict = sm.params.assembly
    # Collect all files belonging to the assembly group
    for sample in assembly_dict.keys():
        for run in assembly_dict[sample]:
            for pair in assembly_dict[sample][run].keys():
                files[pair].append(
                    assembly_dict[sample][run][pair][0])
    # Rename and concatenate reads (required for Metaspades)
    with open(sm.output.R1, 'w') as fh1, open(sm.output.R2, 'w') as fh2, open(
        sm.output.se, 'w') as fhse:
        i = 0
        for f in files["R1"]:
            f2 = files["R2"][i]
            fh1 = rename_records(f, fh1, i)
            fh2 = rename_records(f2, fh2, i)
            i += 1
        for i, f in enumerate(files["se"], start=i):
            fhse = rename_records(f, fhse, i)


def megahit_input(sm):
    """
    Genereate input lists for megahit assembler

    :param sm: snakemake object
    :return:
    """
    files = {"R1": [], "R2": [], "se": []}
    assembly_dict = sm.params.assembly
    for sample in assembly_dict.keys():
        for run in assembly_dict[sample]:
            for pair in assembly_dict[sample][run].keys():
                files[pair].append(assembly_dict[sample][run][pair][0])
    with open(sm.output.R1, 'w') as fh1, \
        open(sm.output.R2, 'w') as fh2, \
        open(sm.output.se, 'w') as fhse:
        fh1.write(",".join(files["R1"]))
        fh2.write(",".join(files["R2"]))
        fhse.write(",".join(files["se"]))


def main(sm):
    if sm.params.assembler == "metaspades":
        metaspades_input(sm)
    elif sm.params.assembler == "megahit":
        megahit_input(sm)


if __name__ == "__main__":
    main(snakemake)
