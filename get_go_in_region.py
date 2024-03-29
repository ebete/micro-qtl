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

chrmap = {
    "1": "NC_015438.1",
    "2": "NC_015439.1",
    "3": "NC_015440.1",
    "4": "NC_015441.1",
    "5": "NC_015442.1",
    "6": "NC_015443.1",
    "7": "NC_015444.1",
    "8": "NC_015445.1",
    "9": "NC_015446.1",
    "10": "NC_015447.1",
    "11": "NC_015448.1",
    "12": "NC_015449.1"
}


def get_peaks(peaks_file):
    """
    Read the QTL peak locations from a CSV file.

    :type peaks_file: str
    :param peaks_file: The location of the CSV file.

    :rtype: list[tuple[str, str, int, int]]
    :return: A list of tuples containing the LOD, chromosome, start, and end.
    """
    lod_peaks = list()
    with open(peaks_file, "r", newline="") as f:
        handle = csv.reader(f, delimiter="\t")
        next(handle)
        for record in handle:
            lod_peaks.append((record[0], record[1], int(record[2]), int(record[3])))
    return lod_peaks


def get_geneid2go(idmapping_file):
    logging.info("Decompressing and loading ID mapping database from %s ...", idmapping_file)
    with lzma.open(idmapping_file, "rb") as f:
        return pickle.load(f)


def load_gff_genes(gff_file):
    logging.info("Parsing genes from %s ...", gff_file)
    gene_records = dict()
    with gzip.open(gff_file, "rt") as gff:
        for line in gff:
            if not line:
                continue
            if line[0] == "#":
                continue
            line = line.strip().split("\t")
            record = {
                "chr": line[0],
                "start": int(line[3]),
                "end": int(line[4]),
                "type": line[2],
                "annotation": line[-1]
            }

            if record["type"] != "gene":
                # skip all non-genes
                continue

            record["annotation"] = {k: v for k, v in (x.split("=") for x in record["annotation"].split(";"))}
            record["xref"] = {k: v for k, v in (x.split(":") for x in record["annotation"]["Dbxref"].split(","))}
            gene_id = record["xref"].get("GeneID")
            if gene_id is not None:
                gene_records[gene_id] = record

    return gene_records


def extract_genes_from_regions(gff_records, chr_id, region_start, region_end):
    logging.info("Extracting genes in %s:%s-%s ...", chr_id, region_start, region_end)

    matching_records = set()
    for gene_id, record in gff_records.items():
        if chr_id is not None and record["chr"] != chr_id:
            # filtering based on chromosome
            continue
        if region_start is not None and record["start"] < region_start:
            # filtering based on region start
            continue
        if region_end is not None and record["end"] > region_end:
            # filtering based on region end
            continue
        matching_records.add(gene_id)

    return matching_records


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
        lookup_table = get_geneid2go(args.mapping_file)
        gff_genes = load_gff_genes(args.gff_file)

        print("lod", "peak_id", "gi", "go", sep="\t")
        peaks = get_peaks(args.peaks_file) #+ [("all", None, None, None)]
        peak_count = 1
        for lod, chromosome, start, end in peaks:
            chromosome = chrmap.get(chromosome, chromosome)
            genes_in_region = extract_genes_from_regions(gff_genes, chromosome, start, end)
            gi_to_go = get_go_terms(genes_in_region, lookup_table)
            for gi, go in gi_to_go.items():
                print(lod, peak_count, gi, "; ".join(go), sep="\t")
            peak_count += 1
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
