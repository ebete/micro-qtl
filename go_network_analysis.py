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

from go_helpers import *

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


def get_go_terms(gi_go_dict):
    all_go = list()
    for go_terms in gi_go_dict.values():
        all_go += go_terms
    return all_go


def propagate_scores(go_tree, go_term, scores, score_to_add):
    for parent in go_tree[go_term].parents:
        scores[parent] = scores.get(parent, 0) + score_to_add
        propagate_scores(go_tree, parent, scores, score_to_add)


def single_network_analysis(go_tree, global_term_occurrence, found_go_terms):
    term_occurrence = get_value_frequency(found_go_terms)
    # relative_term_freq = {k: v / global_term_occurrence[k] for k, v in term_occurrence.items()}

    network_impact = dict()
    for term, occurrence in term_occurrence.items():
        scores = {term: occurrence}
        propagate_scores(go_tree, term, scores, occurrence)
        # if "GO:0008150" not in scores:
            # remove non-BP terms
        # continue
        for lineage_term, propagated_score in scores.items():
            network_impact[lineage_term] = network_impact.get(lineage_term, 0) \
                                           + propagated_score / global_term_occurrence[lineage_term]

    return network_impact


def make_dot_graph(go_tree, term_impact_scores, graph_name="G"):
    digraph_list = list()
    min_v = min(term_impact_scores.values())
    max_v = max(term_impact_scores.values())
    if min_v == max_v:
        max_v += 1

    print(f'digraph "Peak {graph_name}" {{')

    print("rankdir = RL;")
    print("node[shape = ellipse];")
    print("graph[splines = ortho, nodesep = 0.5, bgcolor = white];")
    print()

    for node, impact in term_impact_scores.items():
        traverse_tree(go_tree, digraph_list, node)
        scaled_score = (impact - min_v) / (max_v - min_v)
        hsv_color = f"1.000 1.000 {scaled_score:.3f}"  # if scaled_score > 0.01 else "0.667 1.000 1.000"
        # DOT vertex styling
        print(
            f'"{node}" [style = "filled", fillcolor = "{hsv_color}", fontcolor = "white", '
            f'href="https://www.ebi.ac.uk/QuickGO/term/{node}", target="_blank", '
            f'tooltip="{go_tree[node].go_name}"];'
        )

    print()

    # make DAG
    flattend_tree = set(digraph_list)
    inf_gain = [go_tree[v1].information_content - go_tree[v2].information_content for v1, v2 in flattend_tree]
    max_gain = max(inf_gain)
    min_gain = min(inf_gain)

    for v1, v2 in flattend_tree:
        scaled_inf_gain = (go_tree[v1].information_content - go_tree[v2].information_content - min_gain) / (
                    max_gain - min_gain)
        print(
            f'"{v1}" -> "{v2}" '
            f'[color = "0.800 1.000 {max(0.1, scaled_inf_gain):.3f}"];'
        )

    print("}")


def traverse_tree(go_tree, tuple_list, current_node):
    go_node = go_tree[current_node]

    if not go_node.parents:
        return

    for parent in go_node.parents:
        tuple_list.append((current_node, parent))
        traverse_tree(go_tree, tuple_list, parent)


def show_top(go_tree, term_impact_scores, relative_occurrence, n=10):
    print("go_term", "go_name", "included_in_peaks", "global_fraction_in_peaks", sep="\t")
    for k, v in sorted(term_impact_scores.items(), key=lambda x: (x[1], relative_occurrence[x[0]], x[0]), reverse=True)[
                :n]:
        print(k, go_tree[k].go_name, f"{v:.3f}", f"{relative_occurrence[k]:.3f}", sep="\t")


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

        genome_go_terms = get_go_terms(all_terms["all"])

        global_occurrence = get_value_frequency(genome_go_terms)
        # global_lineage_occurrence = go_lineage_frequencies(go_tree, genome_go_terms)

        lod_occurrence = dict()
        term_occurrence = dict()
        relative_occurrence = dict()
        for lod, gi_to_go in all_terms.items():
            if lod == "all":
                continue
            ancestors = set()
            terms_in_region = [x for x in get_go_terms(gi_to_go) if x in go_tree]
            for term in terms_in_region:
                get_all_ancestors(go_tree, term, ancestors)
                relative_occurrence[term] = relative_occurrence.get(term, 0) + 1 / global_occurrence[term]
            for term in set(terms_in_region):
                term_occurrence[term] = term_occurrence.get(term, 0) + 1 / (len(all_terms) - 1)
            lod_occurrence[lod] = ancestors

            # term_impact = single_network_analysis(go_tree, global_lineage_occurrence, get_go_terms(gi_to_go))
            # if not term_impact:
            #     continue
            #
            # make_dot_graph(go_tree, term_impact, lod)
        show_top(go_tree, term_occurrence, relative_occurrence, n=10)
        # make_dot_graph(go_tree, term_occurrence)
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
