#!/usr/bin/env python

import sys
import pandas as pd
from argparse import ArgumentParser


def add_lower(df, ranks):
    for i in df.index:
        last_known = df.loc[i,ranks[0]]
        for rank in ranks[1:]:
            if df.loc[i,rank] != "Unclassified":
                last_known = df.loc[i,rank]
            else:
                if last_known == "Unclassified":
                    df.loc[i, rank] = last_known
                else:
                    df.loc[i, rank] = "Unclassified.{}".format(last_known)
    return df


def main(args):
    # Keep stats on assignments
    # resolved = cases where sourmash helped resolve assignments
    # transferred = cases where blast-based assignments were overwritten
    # added = cases where assignments from sourmash were added
    # total = total number of contigs
    stats = {'resolved': 0, 'transferred': 0, 'added': 0, 'total': 0}
    df1 = pd.read_csv(args.sourmash, sep=",", header=0, index_col=0)
    stats['total'] = df1.shape[0]
    df2 = pd.read_csv(args.tango, sep="\t", header=0, index_col=0)
    ranks = list(df2.columns)
    ranks.reverse()
    # Only use subset of contigs with matches
    df1 = df1.loc[df1["status"] == "found", df2.columns]
    df1.fillna("Unclassified", inplace=True)

    # Get common set of contigs
    common = set(df1.index).intersection(set(df2.index))
    for contig in common:
        s = df1.loc[contig]
        b = df2.loc[contig]
        for rank in ranks:
            # If sourmash has an assignment at this rank
            if s[rank] != "Unclassified":
                # If blast-based contains 'Unclassified',
                # mark contig as resolved
                if "Unclassified" in b[rank]:
                    stats['resolved'] += 1
                # Otherwise, mark contig as transferred
                else:
                    stats['transferred'] += 1
                # As soon as a contig has been transferred or resolved
                # we can stop the merge
                df2.loc[contig] = df1.loc[contig]
                break
            # If sourmash does not have an assignment at this rank
            else:
                # but blast-based does have an assignment,
                # then the blast-based is more resolved and we can stop
                # trying to merge
                if "Unclassified" not in b[rank]:
                    break
    # Get contigs in sourmash missing from blast
    missing1 = set(df1.index).difference(set(df2.index))
    if len(missing1) > 0:
        stats['added']+=len(missing)
        df2 = pd.concat([df2, df1.loc[missing]])
    df2 = add_lower(df2, df2.columns)
    df2.to_csv(sys.stdout, sep="\t")
    sys.stderr.write("Total:       {}\n".format(stats['total']))
    sys.stderr.write("Resolved:    {}\n".format(stats['resolved']))
    sys.stderr.write("Transferred: {}\n".format(stats["transferred"]))
    sys.stderr.write("Added:       {}\n".format(stats['added']))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("sourmash", type=str,
                        help="Sourmash-based taxonomy assignments")
    parser.add_argument("tango", type=str,
                        help="BLAST-based taxonomy assignments")
    args = parser.parse_args()
    main(args)