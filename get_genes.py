#!/usr/bin/env python3

""""
Copyright (c) 2019 Thom Griffioen
MIT License
"""

import csv
import gzip
import logging
import lzma
import pickle
import sys
import argparse
from pathlib import Path

__description__ = "TBA."
__epilog__ = """
TBA.
"""
__version__ = "2019.9.23"


def get_lod_regions(peaks_file):
    regions = list()
    with open(peaks_file.expanduser().resolve(), "r", newline="") as csv_file:
        handle = csv.reader(csv_file, delimiter='\t')
        next(handle)
        for record in handle:
            regions.append(record)
    return regions


def refseq2chr(mapping_file):
    mapping = dict()
    with open(mapping_file.expanduser().resolve(), "r", newline="") as csv_file:
        handle = csv.reader(csv_file, delimiter='\t')
        next(handle)
        for record in handle:
            mapping[record[0]] = record[1]
    return mapping


def get_genes_from_gff(gff_file, refseq2chr_file, region):
    mapping = refseq2chr(refseq2chr_file)
    matching_records = list()
    with gzip.open(gff_file.expanduser().resolve(), "rt") as gff:
        for line in gff:
            if not line:
                continue
            if line[0] == "#":
                continue
            line = line.strip().split("\t")
            record = dict(
                chr=line[0],
                start=int(line[3]),
                end=int(line[4]),
                type=line[2],
                annotation=line[-1]
            )

            if record["type"] != "gene":
                continue

            find_chr = mapping.get(region[1], region[1])
            if record["chr"] != find_chr or record["start"] < int(region[2]) or record["end"] > int(region[3]):
                continue

            matching_records.append(record)
    return matching_records


def get_geneid2go(idmapping_file):
    logging.info("Decompressing and loading ID mapping database from %s ...", idmapping_file)
    with lzma.open(idmapping_file.expanduser().resolve(), "rb") as f:
        return pickle.load(f)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("idmap_file", metavar="ID_MAP", help="File containing the pickled ID mapping from "
                                                             "make_uniprot_idmapping_db.py")
    parser.add_argument("gff_file", metavar="GENOME_GFF", help="GFF file (gzipped) of the genome containing the genes")
    parser.add_argument("peaks_file", metavar="PEAKS", help="CSV file containing the positions of the QTL's LOD peaks")
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

    regions = get_lod_regions(Path(args.peaks_file))
    geneid2go = get_geneid2go(Path(args.idmap_file))
    logging.debug("\t".join([
        "lod_index",
        "gene_index",
        "chr",
        "gene_id",
        "gene_name",
        "go_terms"
    ]))

    print("lod", "geneid", "length", "go_terms", sep="\t")
    go_in_peak = dict()
    for region in regions:
        logging.debug("region: " + "\t".join([str(x) for x in region]))
        records = get_genes_from_gff(Path(args.gff_file), Path("~/micro_qtl/data/chr2refseq_SL2.40.csv"), region)
        go_per_gene = dict()
        idx = 0
        for record in records:
            anno = record["annotation"].split(";")
            anno_kv = {k: v for k, v in (x.split("=") for x in anno)}
            xref = anno_kv["Dbxref"].split(",")
            xref_kv = {k: v for k, v in (x.split(":") for x in xref)}

            gene_id = xref_kv.get("GeneID", "-")
            go_terms = geneid2go.get(gene_id, list())

            logging.debug("\t".join([
                region[0],
                str(idx),
                record["chr"],
                gene_id,
                anno_kv.get("gene", "-"),
                ";".join(go_terms)
            ]))

            if len(go_terms) > 0:
                go_per_gene[gene_id] = go_per_gene.get(gene_id, list()) + go_terms
                go_in_peak[region[0]] = go_in_peak.get(region[0], list()) + go_terms

            idx += 1
        for gene_id, go_terms in go_per_gene.items():
            print(region[0], gene_id, record["end"]-record["start"], "; ".join(set(go_terms)), sep="\t")

    # for lod_peak, go_terms in go_in_peak.items():
    #     print(lod_peak, "; ".join(set(go_terms)), sep="\t")
