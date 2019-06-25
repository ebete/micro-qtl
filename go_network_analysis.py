#!/usr/bin/env python3

""""
Copyright (c) 2019 Thom Griffioen
MIT License
"""

import argparse
import csv
import logging
import lzma
import pickle
import sys

__description__ = "TBA."
__epilog__ = """
TBA.
"""
__version__ = "2019.6.0"


def read_terms(csv_file):
    lod_gi_to_go = dict()

    logging.info("Reading GI/GO tuples from %s ...", csv_file)
    with open(csv_file, "r", newline="") as f:
        handle = csv.reader(f, delimiter="\t")
        next(handle)
        for record in handle:
            lod = record[0]
            gi_term = record[1]
            go_terms = record[2].split("; ")
            # removes empty splits
            lod_gi_to_go.setdefault(lod, dict())[gi_term] = [x for x in go_terms if len(x) > 0]

    return lod_gi_to_go


def import_go_tree(import_location):
    logging.info("Decompressing and importing GO dictionary from %s ...", import_location)
    with lzma.open(import_location, "rb") as f:
        return pickle.load(f)


def get_go_terms(gi_go_dict):
    all_go = list()
    for go_terms in gi_go_dict.values():
        all_go += go_terms
    return all_go


def get_all_ancestors(go_tree, go_term, ancestors):
    parents = go_tree[go_term].parents
    for parent in parents:
        ancestors.append(parent)
        get_all_ancestors(go_tree, parent, ancestors)


def go_lin_similarity(go_tree, term1, term2):
    """
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.55.1832&rep=rep1&type=pdf

    :param go_tree:
    :param term1:
    :param term2:
    :return: Lin's term similarity score.
    """
    intersecting_ancestors = lowest_common_ancestor(go_tree, term1, term2)
    if not intersecting_ancestors:
        return 0
    lca = intersecting_ancestors[0]
    # get the LCS with the highest IC
    for term in set(intersecting_ancestors):
        if go_tree[lca].information_content < go_tree[term].information_content:
            lca = term

    # calculate Lin's similarity score
    return 2 * go_tree[lca].information_content / \
           (go_tree[term1].information_content + go_tree[term2].information_content)


def lowest_common_ancestor(go_tree, term1, term2):
    """
    Find the lowest common ancestor (LCS) of all paths in the GO DAG.

    :param go_tree:
    :param term1:
    :param term2:
    :return: List of LCS's found on each possible path.
    """
    go_term1 = go_tree[term1]
    go_term2 = go_tree[term2]

    if go_term1 == go_term2:
        return [term1]

    lcs = list()
    # iterate over parents of the most specific node (lower in tree)
    if go_term1.information_content > go_term2.information_content:
        for parent in go_term1.parents:
            subsumer = lowest_common_ancestor(go_tree, parent, go_term2.go_id)
            if not subsumer:
                continue
            lcs += subsumer
    else:
        for parent in go_term2.parents:
            subsumer = lowest_common_ancestor(go_tree, go_term1.go_id, parent)
            if not subsumer:
                continue
            lcs += subsumer

    return lcs


def propagate_scores(go_tree, go_term, scores):
    for child in go_tree[go_term].children:
        # scores[go_term] = scores.get(go_term, 0) + scores.get(child, 0) * (go_tree[child].total_offspring + 1) / \
        # max(1, go_tree[go_term].total_offspring)
        scores[go_term] = scores.get(go_term, 0) + 1 / max(1, go_tree[go_term].total_offspring)
    for parent in go_tree[go_term].parents:
        propagate_scores(go_tree, parent, scores)


def single_network_analysis(go_tree, found_go_terms):
    term_occurrence = dict()
    for term in found_go_terms:
        term_occurrence[term] = term_occurrence.get(term, 0.) + 1

    network_impact = dict()
    for term, occurrence in term_occurrence.items():
        scores = {term: occurrence / max(1, go_tree[term].total_offspring)}
        propagate_scores(go_tree, term, scores)
        if "GO:0008150" not in scores:
            # remove non-BP terms
            continue
        for lineage_term, propagated_score in scores.items():
            network_impact[lineage_term] = network_impact.get(lineage_term, 0) + propagated_score

    return network_impact


def make_dot_graph(go_tree, term_impact_scores):
    digraph_list = list()
    min_v = min(term_impact_scores.values())
    max_v = max(term_impact_scores.values())

    print("digraph G {")

    for node, impact in term_impact_scores.items():
        traverse_tree(go_tree, digraph_list, node)
        # DOT vertex styling
        print('"{}" [style = "filled", fillcolor = "1.000 1.000 {:.3f}" fontcolor = "white"];'
              .format(node, (impact - min_v) / (max_v - min_v)))

    # make DAG
    for v1, v2 in set(digraph_list):
        print('"{}" -> "{}";'.format(v1, v2))

    print("}")


def traverse_tree(go_tree, tuple_list, current_node):
    go_node = go_tree[current_node]

    if not go_node.parents:
        return

    for parent in go_node.parents:
        tuple_list.append((current_node, parent))
        traverse_tree(go_tree, tuple_list, parent)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("go_tree_file", metavar="GO_TREE", help="File containing the pickled GO tree from "
                                                                "create_go_tree.py")
    parser.add_argument("terms_file", metavar="TERMS", help="CSV file containing GI and GO terms of the LOD "
                                                            "peaks and genome")
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
        go_tree = import_go_tree(args.go_tree_file)
        all_terms = read_terms(args.terms_file)
        for lod, gi_to_go in all_terms.items():
            if lod == "all":
                continue
            term_impact = single_network_analysis(go_tree, get_go_terms(gi_to_go))
            make_dot_graph(go_tree, term_impact)
            break
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
