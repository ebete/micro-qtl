#!/usr/bin/env python3

""""
Copyright (c) 2019 Thom Griffioen
MIT License
"""

import argparse
import gzip
import logging
import sys

__description__ = "Create a tree-like structure of the Gene Ontology flat-file."
__epilog__ = """
TBA.
"""
__version__ = "2019.6.0"

go_roots = {
    "biological_process": "GO:0008150",
    "cellular_component": "GO:0005575",
    "molecular_function": "GO:0003674"
}


class GoTerm(object):
    def __init__(self, go_id, go_namespace, child_terms=None):
        self.go_id = go_id
        self.namespace = go_namespace
        self.children = child_terms
        self.total_offspring = 0

    def __hash__(self):
        # we assume that go_id is unique for all instances (should be the case anyway)
        return hash(self.go_id)

    def __repr__(self):
        return 'GoTerm(go_id="{}", go_namespace="{}", total_offspring={})'.format(
            self.go_id, self.namespace, self.total_offspring
        )

    def add_child(self, child_term):
        if not self.children:
            self.children = [child_term]
        else:
            self.children.append(child_term)
        return self


def parse_go_parents(go_file):
    go_map = dict()

    logging.info("Reading GO data file %s ...", go_file)
    with gzip.open(go_file, "rt") as f:
        # fast-forward to first GO term definition
        for line in f:
            if line.startswith("[Term]"):
                break

        read_record = True
        term_data = dict()
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("[Term]"):
                # write term_data to tree
                if term_data.get("is_obsolete", "false") != "true":
                    parent_nodes = [x.split(" ! ")[0] for x in term_data.get("is_a", list())]
                    go_map[term_data["id"]] = parent_nodes
                # initialise parser for next term
                term_data = dict()
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


def traverse_tree(go_parents, lookup_go):
    parents = go_parents.get(lookup_go, list())

    if not parents:
        return

    for parent in parents:
        print('{} -> {}'.format(lookup_go, parent))
        traverse_tree(go_parents, parent)


def reverse_lookup_map(lookup_map):
    reversed_map = dict()
    for child, parents in lookup_map.items():
        for parent in parents:
            reversed_map.setdefault(parent, list()).append(child)
    return reversed_map


def create_go_tree(go_children_map, current_node):
    children = go_children_map.get(current_node)
    new_node = GoTerm(current_node, "biological_process")

    if not children:
        return new_node

    for child in children:
        new_node.add_child(create_go_tree(go_children_map, child))

    return new_node


def calculate_offspring_counts(go_tree):
    if go_tree.children is None:
        return 0

    for child_node in go_tree.children:
        go_tree.total_offspring += 1 + calculate_offspring_counts(child_node)
    return go_tree.total_offspring


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("go_file", metavar="GO", help="Gene Ontology data file (go-basic.obo.gz)")
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
        go_parents = parse_go_parents(args.go_file)
        # for k, v in go_parents.items():
        #     print(k, v, sep="\t")
        go_children = reverse_lookup_map(go_parents)
        go_bp_tree = create_go_tree(go_children, go_roots["biological_process"])
        calculate_offspring_counts(go_bp_tree)
        print(go_bp_tree)
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
