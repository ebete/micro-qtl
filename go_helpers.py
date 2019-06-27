#!/usr/bin/env python3

import pickle
import logging
import lzma

go_roots = {
    "biological_process": "GO:0008150",
    "cellular_component": "GO:0005575",
    "molecular_function": "GO:0003674"
}


class GoTerm(object):
    """
    A doubly linked list of GO records containing some extra metadata of the
    given GO term.
    """

    def __init__(self, go_id, go_name=None, go_def=None):
        """
        Creates a GO term object.

        :type go_id: str
        :param go_id: The unique identifier given to the GO term. This term has
            to be unique and will be used as the hash of this instance.

        :type go_name: str
        :param go_name: The name of the GO term.

        :type go_def: str
        :param go_def: The definition of the GO term.
        """
        if go_id is None:
            raise ValueError("go_id cannot be None.")

        self.go_id = go_id
        self.go_name = go_name
        self.go_def = go_def
        self.children = list()
        self.parents = list()
        self.total_offspring = 0
        self.information_content = 0.

    def __hash__(self):
        # we assume that go_id is unique for all instances (should be the case anyway)
        return hash(self.go_id)

    def __repr__(self):
        return 'GoTerm(go_id="{}", go_name="{}", go_def="{}", children={}, parents={}, total_offspring={}, ' \
               'information_content={})' \
            .format(
            self.go_id, self.go_name, self.go_def, len(self.children), len(self.parents), self.total_offspring,
            self.information_content
        )

    def __str__(self):
        return '{} [{}]'.format(
            self.go_id, self.go_name
        )

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.go_id == other.go_id

    def set_name(self, go_name):
        """
        Set the name of the GO term.

        :type go_name: str
        :param go_name: The new name of the GO term.

        :rtype: GoTerm
        :return: This GoTerm instance.
        """
        self.go_name = go_name
        return self

    def set_definition(self, go_def):
        """
        Set the definition of the GO term.

        :type go_def: str
        :param go_def: The new definition of the GO term.

        :rtype: GoTerm
        :return: This GoTerm instance.
        """
        self.go_def = go_def
        return self

    def add_parent(self, parent_term):
        """
        Add a new parent GO term to this GO term.

        :type parent_term: str
        :param parent_term: The ID of the parent GO term.

        :rtype: GoTerm
        :return: This GoTerm instance.
        """
        if parent_term not in self.parents:
            self.parents.append(parent_term)
        return self

    def add_child(self, child_term):
        """
        Add a new child GO term to this GO term.

        :type child_term: str
        :param child_term: The ID of the child GO term.

        :rtype: GoTerm
        :return: This GoTerm instance.
        """
        if child_term not in self.children:
            self.children.append(child_term)
        return self


def export_go_tree(go_tree, export_location):
    """
    Serialises and compresses the GO tree object into a single file.

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO dictionary to export.

    :type export_location: str
    :param export_location: The location to write the file to.
    """
    logging.info("Compressing and exporting GO dictionary to %s ...", export_location)
    with lzma.open(export_location, "wb") as f:
        pickle.dump(go_tree, f, protocol=pickle.HIGHEST_PROTOCOL)


def import_go_tree(import_location):
    """
    Decompresses and deserialises the given file containing the processed GO
    tree created by create_go_tree.py

    :type import_location: str
    :param import_location: File location of the LZMA compressed and pickled
        object.

    :rtype: dict[str, GoTerm]
    :return: The deserialised object from the file.
    """
    logging.info("Decompressing and importing GO dictionary from %s ...", import_location)
    with lzma.open(import_location, "rb") as f:
        return pickle.load(f)


def go_lin_similarity(go_tree, term1, term2):
    """
    Calculate Lin's similarity score between two GO terms.
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.55.1832&rep=rep1&type=pdf

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO tree dictionary.

    :type term1: str
    :param term1: The first GO term.

    :type term2: str
    :param term2: The second GO term.

    :rtype: float
    :return: Lin's term similarity score.
    """
    intersecting_ancestors = lowest_common_ancestor(go_tree, term1, term2)
    if not intersecting_ancestors:
        return 0.
    lca = intersecting_ancestors.pop()
    # get the LCS with the highest IC
    for term in intersecting_ancestors:
        if go_tree[lca].information_content < go_tree[term].information_content:
            lca = term

    # calculate Lin's similarity score
    return 2 * go_tree[lca].information_content / \
           (go_tree[term1].information_content + go_tree[term2].information_content)


def lowest_common_ancestor(go_tree, term1, term2):
    """
    Find the lowest common ancestor (LCA) of all paths in the GO DAG.

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO tree dictionary.

    :type term1: str
    :param term1: The first GO term.

    :type term2: str
    :param term2: The second GO term.

    :rtype: set[str]
    :return: Set of LCA's found on each possible path.
    """
    go_term1 = go_tree[term1]
    go_term2 = go_tree[term2]

    if go_term1 == go_term2:
        return {term1}

    lca = set()
    # iterate over parents of the most specific node (lower in tree)
    if go_term1.information_content > go_term2.information_content:
        for parent in go_term1.parents:
            subsumer = lowest_common_ancestor(go_tree, parent, go_term2.go_id)
            if not subsumer:
                continue
            lca.update(subsumer)
    else:
        for parent in go_term2.parents:
            subsumer = lowest_common_ancestor(go_tree, go_term1.go_id, parent)
            if not subsumer:
                continue
            lca.update(subsumer)

    return lca


def get_all_ancestors(go_tree, go_term, ancestors):
    """
    Add all ancestor terms of a given GO term to a list.

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO tree dictionary.

    :type go_term: str
    :param go_term: GO term to find all ancestor terms of.

    :type ancestors: set[str]
    :param ancestors: Set where all ancestor terms will be added to. The
        given GO term will also be added to this set.
    """
    ancestors.add(go_term)
    for parent in go_tree[go_term].parents:
        get_all_ancestors(go_tree, parent, ancestors)


def get_value_frequency(values_list):
    """
    Count the frequency of objects in a list.

    :type values_list: list
    :param values_list: List containing the values to count the frequencies of.

    :rtype: dict[_KT, int]
    :return: A dictionary containing the values and their counts.
    """
    occurrences = dict()
    for term in values_list:
        occurrences[term] = occurrences.get(term, 0) + 1
    return occurrences


def go_lineage_frequencies(go_tree, go_terms):
    """
    Calculate how often ancestor terms occur in the given list of GO terms.

    :type go_tree: dict[str, GoTerm]
    :param go_tree: The GO tree dictionary.

    :type go_terms: list[str]
    :param go_terms: List of GO terms to calculate the ancestor frequencies of.

    :rtype: dict[str, int]
    :return: Dictionary containing the frequency of found ancestor GO terms.
    """
    lineage_frequency = dict()
    for term in go_terms:
        if term not in go_tree:
            # GO term has been deprecated
            continue

        term_lineage = set()
        get_all_ancestors(go_tree, term, term_lineage)

        for ancestor in term_lineage:
            lineage_frequency[ancestor] = lineage_frequency.get(ancestor, 0) + 1
    return lineage_frequency
