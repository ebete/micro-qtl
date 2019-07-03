#!/usr/bin/env bash

CORES=${1:-$(grep -c ^processor /proc/cpuinfo)}

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd "${dir}"

snakemake --use-conda -p --cores "${CORES}"
