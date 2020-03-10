localrules: infernal_update_gff, infernal_mask_fasta

rule run_prodigal:
    input:
        opj(config["results_path"],"assembly","{group}","final_contigs.fa")
    output:
        genes=opj(config["results_path"],"annotation","{group}","final_contigs.ffn"),
        faa=opj(config["results_path"],"annotation","{group}","final_contigs.faa"),
        gff=opj(config["results_path"],"annotation","{group}","final_contigs.gff")
    params: "-p meta"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*2
    shell:
        """
        prodigal -i {input} -d {output.genes} -a {output.faa} -o {output.gff} -f gff {params}
        """

rule run_infernal:
    input:
        fastafile=opj(config["results_path"],"assembly","{group}","final_contigs.fa"),
        db=opj(config["infernal_dbpath"],"Rfam.cm"),
        cl=opj(config["infernal_dbpath"],"Rfam.clanin")
    output:
        opj(config["results_path"],"annotation","{group}","final_contigs.rfam")
    threads: 2
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*10
    shell:
        """
        cmscan --cpu {threads} --oskip --rfam --cut_ga --nohmmonly --tblout {output} --fmt 2 \
        --clanin {input.cl} {input.db} {input.fastafile} > /dev/null
        """

def parse_cmout(f):
    with open(f) as fh:
        lines = []
        idnums = {}
        for i, line in enumerate(fh):
            if line[0] == "#": continue
            line = line.rstrip()
            line = line.rsplit()
            indices = [1,2,3,5,9,10,11,14,16,17]
            target_name, target_accession, query, clan, start, end, strand, gc, score, evalue = [line[x] for x in indices]
            try:
                idnum = idnums[query]
            except KeyError:
                idnum = 0
                idnums[query] = idnum
            this_idnum = idnum+1
            idnums[query] = this_idnum
            attributes = ["ID="+query+"ncRNA_"+str(this_idnum),"Name="+target_name,"Accession="+target_accession,"Clan="+clan,"GC="+gc,"E-value="+evalue]
            # seqid, source, type, start, end, score, strand, phase, attributes
            gffline = " ".join([query,"cmscan","ncRNA",start,end,score,strand,".",";".join(attributes)])
            lines.append(gffline)
    return lines

rule infernal_update_gff:
    input:
        gff=opj(config["results_path"],"annotation","{group}","final_contigs.gff"),
        cmout=opj(config["results_path"],"annotation","{group}","final_contigs.rfam"),
    output:
        gff=opj(config["results_path"],"annotation","{group}","final_contigs.infernal.gff")
    run:
        import gffutils
        lines = parse_cmout(input.cmout)
        features = []
        for line in lines:
            feature = gffutils.feature.feature_from_line(line, strict = False, dialect = gffutils.constants.dialect)
            features.append(feature)
        gffdb = gffutils.create_db(input.gff, (input.gff)+".db", force=True)
        gffdb.update(features)
        with open(output.gff, 'w') as fh:
            for feature in gffdb.all_features():
                fh.write(str(feature)+"\n")

rule infernal_mask_fasta:
    input:
        fastafile = opj(config["results_path"],"assembly","{group}","final_contigs.fa"),
        cmout = opj(config["results_path"],"annotation","{group}","final_contigs.rfam")
    output:
        fastafile = opj(config["results_path"],"annotation","{group}","final_contigs.masked.fa")
    run:
        from Bio.SeqIO import parse as seq_parser
        lines = parse_cmout(input.cmout)
        mask = {}
        for line in lines:
            line = line.rsplit()
            indices = [0,3,4,6]
            seqid,start,end,strand = [line[x] for x in indices]
            if strand == "-":
                start,end = end,start
            mask[seqid] = (int(start)-1,int(end))

        with open(output.fastafile, 'w') as fh:
            for record in seq_parser(input.fastafile, "fasta"):
                try:
                    start, end = mask[record.id]
                except KeyError:
                    fh.write(">"+record.description+"\n"+str(record.seq)+"\n")
                    continue
                s = record.seq
                masked = str(s[0:start]+len(s[start:end])*"N"+s[end:])
                fh.write(">"+record.description+"\n"+masked+"\n")
