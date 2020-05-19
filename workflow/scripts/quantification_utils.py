#!/usr/bin/env python


def write_featurefile(sm, score=".", group="gene_id", phase="."):
    with open(sm.input[0], 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for line in fhin:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            items = line.split("\t")
            contig = items[0]
            source = items[1]
            method = items[2]
            start = items[3]
            stop = items[4]
            _strand = items[6]
            geneid = items[8]
            if _strand in ['+1', '1', '+']:
                strand = '+'
            else:
                strand = '-'
            gene_id = "{} {}\n".format(group, geneid.split(";")[0].split("=")[-1])
            fhout.write("\t".join([contig, source, method, start, stop, score,
                                   strand, phase, gene_id]))
    return


def main(sm):
    toolbox = {"write_featurefile": write_featurefile}

    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
