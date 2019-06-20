#!/usr/bin/env python3

""""
Copyright (c) 2019 Thom Griffioen
MIT License
"""

import argparse
import csv
import logging
import sys

__description__ = "TBA."
__epilog__ = """
TBA.
"""
__version__ = "2019.6.0"


def read_terms(csv_file):
    gi_to_go = dict()

    logging.info("Reading GI/GO tuples from %s ...", csv_file)
    with open(csv_file, "r", newline="") as f:
        handle = csv.reader(f, delimiter="\t")
        next(handle)
        for record in handle:
            terms = record[1].split("; ")
            # removes empty splits
            gi_to_go[record[0]] = [x for x in terms if len(x) > 0]

    return gi_to_go


def get_go_terms(gi_go_dict):
    all_go = list()
    for go_terms in gi_go_dict.values():
        all_go += go_terms
    return all_go


def calculate_term_generality(gene_go):
    gi_with_terms = sum([1 for x in gene_go.values() if len(x) > 0])
    go_generality = dict()
    for gi, go_terms in gene_go.items():
        for go in go_terms:
            go_generality[go] = go_generality.get(go, 0) + 1 / gi_with_terms
    return go_generality


def prioritise_genes(region_terms, genome_terms):
    # region_go = get_go_terms(region_terms)
    # genome_go = get_go_terms(genome_terms)
    go_generality = calculate_term_generality(genome_terms)

    prioritised_genes = dict()
    for gi, go_terms in region_terms.items():
        # still have to check if term exists in >50% of QTL regions
        prioritised_terms = [x for x in go_terms if go_generality.get(x, 0) < 0.01]
        prioritised_genes[gi] = prioritised_terms

    return prioritised_genes


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("region_file", metavar="REGION", help="CSV file containing GI and GO terms in the LOD region")
    parser.add_argument("genome_file", metavar="GENOME", help="CSV file containing GI and GO terms in the genome")
    # Optional arguments
    # Standard arguments
    parser.add_argument("-v", "--verbose", help="Increase verbosity level", action="count")
    parser.add_argument("-q", "--silent", help="Suppresses output messages, overriding the --verbose argument",
                        action="store_true")
    parser.add_argument("-l", "--log", help="Set the logging output location",
                        type=argparse.FileType('w'), default=sys.stderr)
    parser.add_argument("-V", "--version", action="version", version=__version__)
    return parser.parse_args()


def set_logging(args):
    log_level = logging.WARNING
    if args.silent:
        log_level = logging.ERROR
    elif not args.verbose:
        pass
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    elif args.verbose == 1:
        log_level = logging.INFO
    logging.basicConfig(format="[%(asctime)s] %(message)s", level=log_level, stream=args.log)
    logging.debug("Setting verbosity level to %s" % logging.getLevelName(log_level))


if __name__ == "__main__":
    args = parse_arguments()
    set_logging(args)

    exitcode = 0
    try:
        region_terms = read_terms(args.region_file)
        genome_terms = read_terms(args.genome_file)
        prio = prioritise_genes(region_terms, genome_terms)
        for k, v in prio.items():
            if len(v) > 0:
                print(k, "; ".join(v), sep="\t")
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
