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
import time

__description__ = "Create a dictionary object for converting between database IDs."
__epilog__ = """
Create a dictionary object for converting between database IDs.
The data file can be downloaded at:

ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz .

Column names:
"""
__version__ = "2019.6.0"

_columns = [
    "uniprotkb_ac",
    "uniprotkb_id",
    "geneid",
    "refseq",
    "gi",
    "pdb",
    "go",
    "uniref100",
    "uniref90",
    "uniref50",
    "uniparc",
    "pir",
    "ncbi_taxon",
    "mim",
    "unigene",
    "pubmed",
    "embl",
    "embl_cds",
    "ensembl",
    "ensembl_trs",
    "ensembl_pro",
    "additional_pubmed"
]
_pretty_columns = [
    "UniProtKB-AC",
    "UniProtKB-ID",
    "GeneID (EntrezGene)",
    "RefSeq",
    "GI",
    "PDB",
    "GO",
    "UniRef100",
    "UniRef90",
    "UniRef50",
    "UniParc",
    "PIR",
    "NCBI-taxon",
    "MIM",
    "UniGene",
    "PubMed",
    "EMBL",
    "EMBL-CDS",
    "Ensembl",
    "Ensembl_TRS",
    "Ensembl_PRO",
    "Additional PubMed"
]

for i in range(len(_columns)):
    __epilog__ += f"{_columns[i]} ({_pretty_columns[i]})\n\n"


def parse_tabfile(uprot_tabfile, from_id, to_id):
    csv.field_size_limit(sys.maxsize)  # needed for that one record (U5Z754) in the idmapping file ...
    decode_db = dict()
    start_time = time.time()
    last_time = time.time()
    status_every = 1e6

    logging.info("Creating a lookup table for %s -> %s",
                 _pretty_columns[_columns.index(from_id)], _pretty_columns[_columns.index(to_id)])

    logging.info("Parsing ~160M records from %s ...", uprot_tabfile)
    with gzip.open(uprot_tabfile, "rt", newline="") as f:
        handle = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE, fieldnames=_columns)
        try:
            for record in handle:
                if len(record) != len(_columns):
                #     logging.error("Invalid number of columns in idmapping file. Offending row %d:\n%s",
                #                   lines_processed, record)
                #     raise Exception("Input idmapping file malformed.")
                    continue

                keys = record[from_id]
                values = record[to_id]
                if not keys or not values:
                    # the key/value is empty for this record
                    pass
                else:
                    keys = keys.split("; ")
                    values = values.split("; ")
                    for k in keys:
                        decode_db[k] = decode_db.get(k, list()) + values

                if handle.line_num % status_every == 0:
                    logging.info("Read %dM ID records in %.1fs (+%.1fs)",
                                 handle.line_num // 1e6, time.time() - start_time, time.time() - last_time)
                    last_time = time.time()
        except:
            logging.debug("Exception triggered while parsing record %d", handle.line_num)
            raise
    return decode_db


def write_db(db_dict, output_file):
    logging.info("Compressing and exporting ID mapping database to %s ...", output_file)
    with lzma.open(output_file, "wb") as f:
        pickle.dump(db_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("idmap_file", metavar="UNIPROT_IDMAP", help="Uniprot idmapping_selected.tab.gz file.")
    parser.add_argument("pickle_output", metavar="DB", help="Output location for the generated ID lookup file.")
    # Optional arguments
    parser.add_argument("-f", "--from", metavar="ID", dest="from_id", help="The key of the lookup table.",
                        default="uniprotkb_ac")
    parser.add_argument("-t", "--to", metavar="ID", dest="to_id", help="The value of the lookup table.", default="go")
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
        db = parse_tabfile(args.idmap_file, args.from_id, args.to_id)
        write_db(db, args.pickle_output)
    # except Exception as ex:
    #     exitcode = 1
    #     logging.error(ex)
    #     logging.debug(format_exc())
    finally:
        logging.debug("Shutting down logging system ...")
        logging.shutdown()
    sys.exit(exitcode)
