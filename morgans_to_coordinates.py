#!/usr/bin/env python3

import argparse
import csv
import logging
import sys

from pathlib import Path

__description__ = "TBA."
__epilog__ = """
TBA.
"""
__version__ = "2019.8.0"

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


def read_gen(fname):
    # I f*cking hate this code ...
    with open(fname, "r") as genfile:
        snps = next(genfile)
        snps = [x.split("-")[0] for x in snps.strip().split(",")[1:]]
        chrom = next(genfile)
        chrom = chrom.strip().split(",")[1:]
        cm = next(genfile)
        cm = list(map(float, cm.strip().split(",")[1:]))
        next(genfile)

        rils = dict()
        for ril in genfile:
            ril = ril.strip().split(",")
            rils[ril[0]] = ril[1:]

        rm_idx = [i for i in range(len(snps)) if "seq" in snps[i]]
        for i in rm_idx[::-1]:
            del chrom[i]
            del snps[i]
            del cm[i]
            for ril in rils:
                del rils[ril][i]

        d = dict()
        for i in range(len(snps)):
            d.setdefault(chrom[i], list()).append([snps[i], cm[i]])

        return d


def read_markers(fname):
    markers = list()
    with open(fname, "r", newline="") as peaksfile:
        handle = csv.reader(peaksfile, delimiter="\t")
        next(handle)
        for record in handle:
            markers.append(record)
    return markers


def parse_arguments():
    parser = argparse.ArgumentParser(description=__description__, epilog=__epilog__)
    # Required arguments
    parser.add_argument("marker_file", metavar="GENFILE", help="File containing the genetic marker map")
    parser.add_argument("peaks_file", metavar="PEAKS", help="File containing the significant markers and CI")
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

    marker_file = Path(args.marker_file)
    peaks_file = Path(args.peaks_file)
    # peaks_name = peaks_file.name.partition(".")[0]

    d = read_gen(marker_file)
    markers = read_markers(peaks_file)

    #print("chr", "cM", "pos_low", "pos", "pos_hi", sep="\t")
    print("lod", "chr", "start", "end", sep="\t")
    for marker in markers:
        snps = d[marker[1]]
        for idx in range(len(snps)):
            snp_prev = snps[max(0, idx-1)]
            snp = snps[idx]
            snp_next = snps[min(idx+1, len(snps)-1)]
            if snp[1] < float(marker[2]):
                continue
            #print(marker[0], marker[1], f"chr{marker[0]}:{snp_prev[0]}", f"chr{marker[0]}:{snp[0]}", f"chr{marker[0]}:{snp_next[0]}", sep="\t")
            print(marker[0], chrmap[marker[1]], snp_prev[0], snp_next[0], sep="\t")
            break
