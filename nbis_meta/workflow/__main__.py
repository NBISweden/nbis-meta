#!/usr/bin/env python

from argparse import ArgumentParser
from snakemake import snakemake
from importlib import resources


def main():
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", type=str,
                        default="config/config.yaml", help="Configuration file")
    parser.add_argument("-n", "--dryrun", action="store_true")
    args, sm_args = parser.parse_known_args()
    with resources.path("workflow", "Snakefile") as snakefile:
        s = snakemake(snakefile, dryrun=args.dryrun, configfiles=[args.config])


if __name__ == "__main__":
    main()
