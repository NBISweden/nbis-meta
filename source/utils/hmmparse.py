#!/usr/bin/env python

import sys, pandas as pd, logging
from argparse import ArgumentParser

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def calculate_overlap(ali1, ali2):
    if len(ali1) < len(ali2):
        shortest = len(ali1)
    else:
        shortest = len(ali2)
    # 'overlap' is the number of positions that overlap
    overlap = set(ali1).intersection(set(ali2))
    # divide overlap by the shortest alignment to get overlap fraction, if greater than the specified fraction \
    # then the alignment under investigation is overlapping a hit with a higher score and we should not store it
    overlap_frac = float(len(overlap)) / shortest
    return overlap_frac


def test_calculate_overlap():
    # ---------- #
    #      ----- #
    ali1 = range(1,101)
    ali2 = range(51,101)
    assert calculate_overlap(ali1,ali2) == 1.0
    # -------- #
    #        -- #
    ali3 = range(91,111)
    assert calculate_overlap(ali1, ali3) == 0.5
    # -------- #
    #           ---- #
    ali4 = range(121,201)
    assert calculate_overlap(ali1, ali4) == 0.0


def checkoverlap(r, frac=0.1):
    # alis is a list of ranges defined by the 'from' and 'to' columns of r
    alis = []
    # store is a list of indexes in r
    store = []
    # iterate through the lines in r
    for i in range(0, len(r)):
        # get the 'from' (f) and 'to' (t) coordinates
        f = int(r.iloc[i, 6])
        t = int(r.iloc[i, 7])
        # generate a list of positions between from and to
        this_ali = range(f, t + 1)
        this_ali_len = len(this_ali)
        # If this is the first line, store it (assume lines are sorted by domain score)
        if len(alis) == 0:
            store.append(i)
            alis.append(this_ali)
            continue
        # For each alignment that has been stored, check the amount of overlap with the alignment under investigation
        overlapping = False
        for ali in alis:
            overlap_frac = calculate_overlap(ali, this_ali)
            if overlap_frac > frac:
                overlapping = True
                break
        if not overlapping:
            alis.append(this_ali)
            store.append(i)
    return store


def check_output_type(line):
    header = line.rstrip()
    header_len = len(header.rsplit())
    if header_len == 23:
        return 'domtblout'
    elif header_len == 19:
        return 'tblout'
    else:
        sys.exit("Header length {i} doesn't match either --tblout or --domtblout format. Unknown file input".format(i = header_len))


def readlines(fh):
    r = []
    output_type = "tblout"
    for i, line in enumerate(fh):
        if line[0] == "#" and i == 2:
            output_type = check_output_type(line)
        if line[0] == "#":
            continue
        line = line.rstrip()
        line = line.rsplit()
        if output_type == "domtblout":
            # target_name, query_accession, e-value, score, i-evalue, domain_score, from, to
            columns = ["target_name", "query_accession", "evalue", "score", "ievalue", "domain_score", "from", "to"]
            items = [line[0], line[4], float(line[6]), float(line[7]), float(line[12]), float(line[13]), int(line[17]), int(line[18])]
            if items[1] == "-":
                items[1] = line[3]
        elif output_type == "tblout":
            # target_name, query_accession, e-value, score
            columns = ["target_name", "query_accession", "evalue", "score"]
            items = [line[0],line[3],float(line[4]),float(line[5])]
            if items[1] == "-":
                items[1] = line[2]
        r.append(items)
    df = pd.DataFrame(r)
    df.columns = columns
    return df, output_type


def parse_multi(df,args):

    # If an "independent E-value" was specified, filter domains by this e-value before counting
    if args.ievalue:
        logging.info("Filtering domains with independent e-value > {e}".format(e=args.ievalue))
        l = len(df)
        df = df.loc[df.ievalue<args.ievalue]
        logging.info("Removed {i} hits".format(i=l - len(df)))

    # Count number of potential hits per protein
    hit_counts = pd.DataFrame(df.groupby("target_name").count().iloc[:, 0])
    hit_counts.columns = ["hits"]

    # Get single hit proteins
    singles = hit_counts.loc[hit_counts.hits==1].index
    ann = df.loc[df.target_name.isin(singles)]

    # Get potential multi hit proteins
    multis = hit_counts.loc[hit_counts.hits>1].index

    logging.info("Parsing {r} records for {m} potential multi-hit targets".format(r=len(df.loc[df.target_name.isin(multis)]), m=len(multis)))
    # Iterate targets with potential multi hits and store highest scoring domains that do not overlap more than
    # the specified threshold
    df_domsort = df.sort_values("domain_score", ascending=False)
    total_targets = len(multis)
    for t, target in enumerate(multis, start=1):
        if t % 1000 == 0:
            logging.info("Parsed {t} of {T}".format(t=t,T=total_targets))
        target_hits_sort = df_domsort.loc[df_domsort.target_name==target]
        hits_to_store = checkoverlap(target_hits_sort, args.overlap)
        ann = pd.concat([ann, target_hits_sort.iloc[hits_to_store]])

    # Sort by start position on protein
    ann.sort_values(["target_name", "from"], ascending=True, inplace=True)
    return ann


def file_to_df(f):
    with open(f) as fh:
        df, output_type = readlines(fh)
    return df, output_type


def main():
    parser = ArgumentParser('''Simple parser that selects the highest scoring HMM per query. \
            Users can set an e-value cutoff for filtering.''')
    parser.add_argument("-i", "--infile", required=True,
            help="hmmer tabular output file")
    parser.add_argument("-e", "--evalue", type=float, required=False,
            help="e-value cutoff to use")
    parser.add_argument("--ievalue", type=float,
            help="Use this e-value to filter domains when parsing multiple hits per proteins. This is a stringent \
            measure of how reliable particular domains may be.")
    parser.add_argument("-o", "--outfile",
            help="Write parsed result to outfile. Defaults to stdout.")
    parser.add_argument("--multi", action="store_true",
            help="Store multiple hits per query. Set domain overlap fraction using '-O'")
    parser.add_argument("-O", "--overlap", type=float, default=0.1,
            help="Maximum allowed overlap (as fraction of shortest alignment compared) between two domain hits. \
            Default is 0.1 (10%). Ignored unless run with '--multi'")

    args = parser.parse_args()

    multi = args.multi

    logging.info("Parsing input file {f}".format(f=args.infile))
    df, output_type = file_to_df(args.infile)

    if args.multi and output_type == "tblout":
        logging.info("WARNING: --multi specified but results appear to be in 'tblout' format. Only best hits will be parsed")
        multi = False

    # Filter by e-value
    if args.evalue:
        logging.info("Filtering hits with e-value > {e}".format(e=args.evalue))
        l = len(df)
        df = df.loc[df.evalue<args.evalue]
        logging.info("Removed {i} hits".format(i=l-len(df)))

    # Sort by score
    df.sort_values("score",ascending=False,inplace=True)

    # If not parsing multiple HMM hits per sequence, just select first best hit
    if multi:
        logging.info("Parsing multiple domains per protein")
        ann = parse_multi(df, args)
        ann = ann.loc[:, ["target_name", "query_accession", "evalue", "score", "ievalue", "domain_score", "from", "to"]]
    else:
        logging.info("Parsing {r} records into single best hits per protein".format(r=len(df)))
        ann = df.groupby("target_name").first().reset_index()
        ann = ann.loc[:,["target_name", "query_accession", "evalue", "score"]]

    # Write results
    if not args.outfile:
        fh = sys.stdout
    else:
        fh = open(args.outfile, 'w')
    ann.to_csv(fh,sep="\t", index=False)

if __name__ == '__main__':
    main()
