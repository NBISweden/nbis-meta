#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys


def main():
    p = ArgumentParser()
    p.add_argument("-i", "--infile", required=True,
                   help="Infile table")
    p.add_argument("-c", "--column", default=6,
                   help="Column to sort unique (1-based)")
    p.add_argument("-s", "--sort", default=15,
                   help="Sort table by column first (1-based)")
    args = p.parse_args()

    column = args.column - 1

    df = pd.read_table(args.infile, header=0, skiprows=1, dtype=str)
    if args.sort:
        sort_column = args.sort - 1
        df.sort_values(df.columns[sort_column], ascending = False, inplace = True)

    df = df.groupby(df.columns[column]).first().reset_index()

    df.to_csv(sys.stdout, sep="\t", index=False, header=True)


if __name__ == '__main__':
    main()
