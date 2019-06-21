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
    def __init__(self, go_id, go_name=None, go_def=None):
        if go_id is None:
            raise ValueError("go_id cannot be None.")

        self.go_id = go_id
        self.go_name = go_name
        self.go_def = go_def
        self.children = list()
        self.parents = list()
        self.total_offspring = 0

    def __hash__(self):
        # we assume that go_id is unique for all instances (should be the case anyway)
        return hash(self.go_id)

    def __repr__(self):
        return 'GoTerm(go_id="{}", go_name="{}", go_def="{}", children={}, parents={}, total_offspring={})'.format(
            self.go_id, self.go_name, self.go_def, len(self.children), len(self.parents), self.total_offspring
        )

    def __str__(self):
        return '{} [{}]'.format(
            self.go_id, self.go_name
        )

    def set_name(self, go_name):
        self.go_name = go_name
        return self

    def set_definition(self, go_def):
        self.go_def = go_def
        return self

    def add_parent(self, parent_term):
        if parent_term not in self.parents:
            self.parents.append(parent_term)
        return self

    def add_child(self, child_term):
        if child_term not in self.children:
            self.children.append(child_term)
        return self


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
    term_data = go_tree[go_term]

    # fix against calculating offspring multiple times
    if term_data.total_offspring == 0:
        for child in term_data.children:
            term_data.total_offspring += offspring_calculation(go_tree, child) + 1

    return term_data.total_offspring


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
        go_tree = parse_go_tree(args.go_file)
        # make_plot(go_tree, "GO:2001317", "/dev/null")
        for k, v in go_roots.items():
            offspring_calculation(go_tree, v)
            logging.info("Total nodes under %s [%s]: %d", v, k, go_tree[v].total_offspring)
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
