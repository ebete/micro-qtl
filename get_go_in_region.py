#!/usr/bin/env python3

""""
Copyright (c) 2019 Thom Griffioen
MIT License
"""

import argparse
import csv
import gzip
import logging
import lzma
import pickle
import sys

__description__ = "TBA."
__epilog__ = """
TBA.
"""
__version__ = "2019.6.0"


def get_peaks(peaks_file):
    lod_peaks = list()
    with open(peaks_file, "r", newline="") as f:
        handle = csv.reader(f, delimiter="\t")
        next(handle)
        for record in handle:
            lod_peaks.append(tuple(record))
    return lod_peaks


def get_geneid2go(idmapping_file):
    logging.info("Decompressing and loading ID mapping database from %s ...", idmapping_file)
    with lzma.open(idmapping_file, "rb") as f:
        return pickle.load(f)


def extract_genes_from_regions(gff_file, chr_id, region_start, region_end):
    logging.info("Extracting genes in %s:%s-%s from %s ...", chr_id, region_start, region_end, gff_file)
    matching_records = list()

    with gzip.open(gff_file, "rt") as gff:
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
                # skip all non-genes
                continue
            if chr_id is not None and record["chr"] != chr_id:
                # filtering based on chromosome
                continue
            if region_start is not None and record["start"] < region_start:
                # filtering based on region start
                continue
            if region_end is not None and record["end"] > region_end:
                # filtering based on region end
                continue

            record["annotation"] = {k: v for k, v in (x.split("=") for x in record["annotation"].split(";"))}
            record["xref"] = {k: v for k, v in (x.split(":") for x in record["annotation"]["Dbxref"].split(","))}

            matching_records.append(record)
    return set(x["xref"].get("GeneID") for x in matching_records)


def get_go_terms(gene_ids, lookup_table):
    return {gi: lookup_table.get(gi, list()) for gi in gene_ids}


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("peaks_file", metavar="PEAKS", help="CSV file containing the positions of the LOD peaks of "
                                                            "the QTL")
    # Optional arguments
    # parser.add_argument("-c", "--chromosome", metavar="CHR", help="The chromosome ID of the region", default=None)
    # parser.add_argument("-s", "--start", metavar="START", help="The start position of the region", type=int,
    #                     default=None)
    # parser.add_argument("-e", "--end", metavar="END", help="The end position of the region", type=int, default=None)
    parser.add_argument("-g", "--gff", metavar="FILE", dest="gff_file",
                        help="GFF file (gzipped) of the genome containing the genes", default="genome.gff.gz")
    parser.add_argument("-m", "--map", metavar="FILE", dest="mapping_file",
                        help="ID mapping file generated by make_uniprot_idmapping_db.py", default="mapping.pickle.xz")
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
        print("lod", "gi", "go", sep="\t")
        lookup_table = get_geneid2go(args.mapping_file)
        peaks = get_peaks(args.peaks_file) + [("all", None, None, None)]
        for lod, chromosome, start, end in peaks:
            if lod == "all":
                genes_in_region = extract_genes_from_regions(args.gff_file, None, None, None)
            else:
                genes_in_region = extract_genes_from_regions(args.gff_file, chromosome, int(start), int(end))
            gi_to_go = get_go_terms(genes_in_region, lookup_table)
            for gi, go in gi_to_go.items():
                print(lod, gi, "; ".join(go), sep="\t")
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
