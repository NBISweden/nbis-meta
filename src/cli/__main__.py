#!/usr/bin/env python

from argparse import ArgumentParser
import snakemake


def main():
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", type=str,
                        default="config/config.yaml",
                        help="Config file")
    parser.add_argument("-n", "--dryrun", action="store_true")
    args, sm_args = parser.parse_known_args()
    snakemake.snakemake("workflow/Snakefile", dryrun=args.dryrun,
                        configfiles=[args.config])


if __name__ == "__main__":
    main()
