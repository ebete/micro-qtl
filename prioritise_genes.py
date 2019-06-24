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


def import_go_tree(import_location):
    logging.info("Decompressing and importing GO dictionary from %s ...", import_location)
    with lzma.open(import_location, "rb") as f:
        return pickle.load(f)


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


def prioritise_genes(region_terms, genome_terms, generality_cutoff=0.01):
    # region_go = get_go_terms(region_terms)
    # genome_go = get_go_terms(genome_terms)
    go_generality = calculate_term_generality(genome_terms)

    prioritised_genes = dict()
    for gi, go_terms in region_terms.items():
        # still have to check if term exists in >50% of QTL regions
        prioritised_terms = [x for x in go_terms if go_generality.get(x, 0) < generality_cutoff]
        prioritised_genes[gi] = prioritised_terms

    return prioritised_genes


def prioritise_goterms(region_terms, genome_terms, generality_cutoff=0.01):
    prio = prioritise_genes(region_terms, genome_terms, generality_cutoff)
    go_terms = list()
    for k, v in prio.items():
        go_terms += v
    return set(go_terms)


def term_overrepresentation(go_tree, genome_terms, region_terms):
    genome_counts = go_lineage_occurrence(go_tree, genome_terms)
    region_counts = go_lineage_occurrence(go_tree, region_terms)
    overrepresented_terms = prioritise_goterms(region_terms, genome_terms)

    informative_terms = list()
    sorted_counts = sorted(region_counts.items(), key=lambda x: x[1], reverse=True)  # sort tuple by occurrence
    for go, count in sorted_counts:
        if go not in overrepresented_terms:
            continue

        term_ic = go_tree[go].information_content
        genome_fraction = count / genome_counts[go]

        if term_ic < 6.64:
            continue

        print(go, count, genome_counts[go], genome_fraction, "{:.2f}".format(term_ic), sep="\t")
        informative_terms.append(go)

    # for g1, g2 in combinations(informative_terms, 2):
    #     lin_sim = go_lin_similarity(go_tree, g1, g2)
    #     if lin_sim > 0.9:
    #         print(g1, g2, lin_sim, sep="\t")
    return informative_terms


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


def go_lineage_occurrence(go_tree, terms):
    term_counts = {}
    for gterm in get_go_terms(terms):
        if gterm not in go_tree:
            continue
        term_counts[gterm] = term_counts.get(gterm, 0) + 1
        term_ancestors = list()
        get_all_ancestors(go_tree, gterm, term_ancestors)
        for ancestor in term_ancestors:
            term_counts[ancestor] = term_counts.get(ancestor, 0) + 1
    return term_counts


def get_all_ancestors(go_tree, go_term, ancestors):
    parents = go_tree[go_term].parents
    for parent in parents:
        ancestors.append(parent)
        get_all_ancestors(go_tree, parent, ancestors)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("go_tree_file", metavar="FILE", help="File containing the pickled GO tree from "
                                                             "create_go_tree.py")
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
        go_tree = import_go_tree(args.go_tree_file)

        region_terms = read_terms(args.region_file)
        genome_terms = read_terms(args.genome_file)

        informative = term_overrepresentation(go_tree, genome_terms, region_terms)

        # prio = prioritise_genes(region_terms, genome_terms)
        # for k, v in prio.items():
        #     if len(v) > 0:
        #         print(k, "; ".join(v), sep="\t")
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
