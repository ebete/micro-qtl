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
from math import log2

from get_go_in_region import get_peaks
from go_helpers import *

__description__ = "TBA."
__epilog__ = """
TBA.
"""
__version__ = "2019.7.0"


def run_analysis(args):
    # process input files
    lod_peaks = {x[0]: x[1:] for x in get_peaks(args.peaks_file)}
    go_tree = import_go_tree(args.go_tree_file)
    all_terms = read_terms(args.terms_file)

    # split genome terms from peak terms
    genome_go_terms = dict()
    peak_go_terms = list()
    terms_in_peak = dict()
    total_peak_nt_length = 0
    for lod, gi_to_go in all_terms.items():
        if lod == "all":  # terms of the entire genome
            genome_go_terms = get_go_terms(gi_to_go)
            continue

        terms_in_peak[lod] = get_go_terms(gi_to_go)
        peak_go_terms += get_go_terms(gi_to_go)
        total_peak_nt_length += lod_peaks[lod][2] - lod_peaks[lod][1]

    # calculate background term rates
    genome_nt_length = 800e6
    global_occurrence = get_value_frequency(genome_go_terms)
    global_density = {k: v / genome_nt_length for k, v in global_occurrence.items()}
    global_lineage_occurrence = go_lineage_frequencies(go_tree, genome_go_terms)
    global_lineage_density = {k: v / genome_nt_length for k, v in global_lineage_occurrence.items()}
    # peak region term rates
    local_occurrence = get_value_frequency(peak_go_terms)
    local_density = {k: v / total_peak_nt_length for k, v in local_occurrence.items()}
    local_lineage_occurrence = go_lineage_frequencies(go_tree, peak_go_terms)
    local_lineage_density = {k: v / total_peak_nt_length for k, v in local_lineage_occurrence.items()}

    # compare background rates to QTL rates
    lfd_density = {k: log2(v / global_density[k]) for k, v in local_density.items()}
    lfd_lineage_density = {k: log2(v / global_lineage_density[k]) for k, v in local_lineage_density.items()}

    # show_top(go_tree, lfd_density, n=10, descending=False)
    show_top(go_tree, lfd_lineage_density, n=10, descending=True)

    pass


def read_terms(csv_file):
    """
    Read the GI terms and their associated GO terms from the given CSV file.

    :type csv_file: str
    :param csv_file: The CSV file containing the GI and GO terms.

    :rtype: dict[str, dict[str, list[str]]]
    :return: A dictionary containing a list of all the GO terms of a GI per
        found LOD peak.
    """
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
    """
    Extract all GO terms (with duplicates) from the given GI to GO term
    dictionary.

    :type gi_go_dict: dict[str, list[str]]
    :param gi_go_dict: Dictionary containing the GO terms associated to the GI.

    :rtype: list[str]
    :return: A list containing all found GO terms.
    """
    all_go = list()
    for go_terms in gi_go_dict.values():
        all_go += go_terms
    return all_go


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
    """
    Go up the GO tree starting at the given node while adding all connections
    to the given list.

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO tree dictionary.

    :type tuple_list: list[tuple]
    :param tuple_list: List containing tuples of all found linked terms
        (child -> parent).

    :type current_node: str
    :param current_node: The node to get all the parent terms from.
    """
    go_node = go_tree[current_node]

    if not go_node.parents:
        return

    for parent in go_node.parents:
        tuple_list.append((current_node, parent))
        traverse_tree(go_tree, tuple_list, parent)


def show_top(go_tree, term_scores, n=None, descending=True):
    """
    Print the top n scoring terms in a tab-delimited table.

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO tree dictionary.

    :type term_scores: dict[str, float]
    :param term_scores: Dictionary containing the GO terms and the associated
        score.

    :type n: int
    :param n: Show only the top n terms.

    :type descending: bool
    :param descending: Print the GO term scores in descending order.
    """
    print("go_term", "go_name", "scores", sep="\t")
    for k, v in sorted(term_scores.items(), key=lambda x: (x[1], x[0]), reverse=descending)[:n]:
        print(k, go_tree[k].go_name, f"{v:.3f}", sep="\t")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("go_tree_file", metavar="GO_TREE", help="File containing the pickled GO tree from "
                                                                "create_go_tree.py")
    parser.add_argument("terms_file", metavar="TERMS", help="CSV file containing GI and GO terms of the LOD "
                                                            "peaks and genome")
    parser.add_argument("peaks_file", metavar="PEAKS", help="CSV file containing the positions of the LOD peaks of "
                                                            "the QTL")
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
        run_analysis(args)
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
