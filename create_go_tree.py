#!/usr/bin/env python3

""""
Copyright (c) 2019 Thom Griffioen
MIT License
"""

import argparse
import gzip
import logging
import lzma
import pickle
import sys
from math import log2

from go_helpers import *

__description__ = "Create a tree-like structure of the Gene Ontology flat-file."
__epilog__ = """
TBA.
"""
__version__ = "2019.6.0"


def parse_go_tree(go_file):
    go_map = dict()

    logging.info("Reading GO data file %s ...", go_file)
    with gzip.open(go_file, "rt") as f:
        # fast-forward to first GO term definition
        while not next(f).startswith("[Term]"):
            pass

        read_record = True
        term_data = dict()
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("[") and read_record:
                # write term_data to dictionary
                if term_data.get("is_obsolete", "false") != "true":
                    parent_nodes = [x.split(" ! ")[0] for x in term_data.get("is_a", list())]
                    node_id = term_data["id"]

                    go_map.setdefault(node_id, GoTerm(node_id)) \
                        .set_name(term_data.get("name")) \
                        .set_definition(term_data.get("def"))

                    for parent in parent_nodes:
                        go_map[node_id].add_parent(parent)
                        go_map.setdefault(parent, GoTerm(parent)).add_child(node_id)
                # initialise parser for next term
                term_data = dict()

            if line.startswith("[Term]"):
                read_record = True
                continue
            elif line.startswith("["):
                # skip records that do not define GO term
                read_record = False
                continue

            if not read_record:
                continue

            key, value = line.split(": ", maxsplit=1)
            if key in ["synonym", "is_a", "consider", "subset", "alt_id", "xref"]:
                # key may occur multiple times
                term_data.setdefault(key, list()).append(value)
            else:
                term_data[key] = value

    return go_map


def traverse_tree(go_tree, lookup_go):
    go_node = go_tree[lookup_go]

    if not go_node.parents:
        return

    for parent in go_node.parents:
        print('"{}" -> "{}";'.format(go_node, go_tree[parent]))
        traverse_tree(go_tree, parent)


def make_plot(go_tree, lookup_go, output_fig):
    print("digraph {")
    traverse_tree(go_tree, "GO:2001317")
    print("}")


def offspring_calculation(go_tree, go_term):
    """
    Calculate the total number of nodes under a term by using a recursive
    depth-first tree traversal approach.

    :type go_tree: dict
    :param go_tree: The GO table containing all terms.

    :type go_term: str
    :param go_term: The GO node to calculate the offspring count of.

    :rtype: int
    :return: The number of offspring of the given node.
    """
    term_data = go_tree[go_term]

    # fix against calculating offspring multiple times
    if term_data.total_offspring == 0:
        for child in term_data.children:
            term_data.total_offspring += offspring_calculation(go_tree, child) + 1

    return term_data.total_offspring


def information_content_calculation(go_tree, go_term, total_terms):
    term_data = go_tree[go_term]
    term_data.information_content = -log2((term_data.total_offspring + 1) / total_terms)
    for child in term_data.children:
        information_content_calculation(go_tree, child, total_terms)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("go_file", metavar="GO", help="Gene Ontology data file (go-basic.obo.gz)")
    # Optional arguments
    parser.add_argument("-o", "--output", metavar="FILE", dest="tree_out", help="Output location of the GO tree")
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
        go_tree = parse_go_tree(args.go_file)
        # make_plot(go_tree, "GO:2001317", "/dev/null")

        for k, v in go_roots.items():
            offspring_calculation(go_tree, v)
            information_content_calculation(go_tree, v, go_tree[v].total_offspring + 1)
            logging.info("Total nodes under %s [%s]: %d", v, k, go_tree[v].total_offspring)

        if args.tree_out is not None:
            export_go_tree(go_tree, args.tree_out)
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
