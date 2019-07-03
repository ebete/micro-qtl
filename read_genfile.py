#!/usr/bin/env python3

import sys
import random

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

        return chrom, snps, cm, rils


def get_random_interval(chrom, snps, cm, rils, interval_cm):
    rand_snp = random.randrange(len(snps))
    rand_chr = chrom[rand_snp]
    cm1 = cm[rand_snp]
    while cm[len(chrom) - chrom[::-1].index(rand_chr) - 1] - cm1 < interval_cm:
        rand_snp = random.randrange(len(snps))
        rand_chr = chrom[rand_snp]
        cm1 = cm[rand_snp]

    offset = 0
    cm2 = cm[rand_snp + offset]
    while abs(cm1 - cm2) < interval_cm:
        offset += 1
        cm2 = cm[rand_snp + offset]

    return rand_snp, rand_snp + offset


if __name__ == "__main__":
    chrom, snps, cm, rils = read_gen(sys.argv[1])

    print("lod_id", "chr", "start", "end", sep="\t")
    for i in range(1, int(sys.argv[2]) + 1):
        snp1, snp2 = get_random_interval(chrom, snps, cm, rils, 10.0)
        print(i, chrmap[chrom[snp1]], snps[snp1], snps[snp2], sep="\t")
