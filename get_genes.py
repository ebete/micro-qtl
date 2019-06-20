#!/usr/bin/env python3

import csv
import gzip
import lzma
import pickle
import logging


# from Bio import Entrez


# Entrez.api_key = "f33ee6aed439ad8d9382a53862e70b02cc08"
# Entrez.email = "t.griffioen@nioo.knaw.nl"


def get_lod_regions(peaks_file):
    regions = list()
    with open(peaks_file, "r", newline="") as csv_file:
        handle = csv.reader(csv_file, delimiter='\t')
        next(handle)
        for record in handle:
            regions.append(record)
    return regions


def refseq2chr(mapping_file):
    mapping = dict()
    with open(mapping_file, "r", newline="") as csv_file:
        handle = csv.reader(csv_file, delimiter='\t')
        next(handle)
        for record in handle:
            mapping[record[0]] = record[1]
    return mapping


def get_genes_from_gff(gff_file, refseq2chr_file, region):
    mapping = refseq2chr(refseq2chr_file)
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
                continue

            find_chr = mapping.get(region[1], "-")
            if record["chr"] != find_chr or record["start"] < int(region[2]) or record["end"] > int(region[3]):
                continue

            matching_records.append(record)
    return matching_records


def get_geneid2go(idmapping_file):
    logging.info("Decompressing and loading ID mapping database from %s ...", idmapping_file)
    with lzma.open(idmapping_file, "rb") as f:
        return pickle.load(f)


if __name__ == "__main__":
    logging.basicConfig(format="[%(asctime)s] %(message)s", level=logging.DEBUG)

    regions = get_lod_regions("ril_lod_peaks.csv")
    geneid2go = get_geneid2go("../db/geneid2go.pickle.xz")
    logging.debug("\t".join([
        "lod_index",
        "gene_index",
        "chr",
        "gene_id",
        "gene_name",
        "go_terms"
    ]))

    go_in_peak = dict()
    for region in regions:
        records = get_genes_from_gff("../db/GCF_000188115.4_SL3.0_genomic.gff.gz", "refseq2chr.csv", region)
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
            print("{}\t{}".format(gene_id, "; ".join(set(go_terms))))

    for lod_peak, go_terms in go_in_peak.items():
        print("Peak {}\t{}".format(lod_peak, "; ".join(set(go_terms))))
