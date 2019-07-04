#!/usr/bin/env bash

CORES=${1:-$(grep -c ^processor /proc/cpuinfo)}

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd "${dir}"

check_updates() {
	UPSTREAM=${1:-'@{u}'}
	LOCAL=$(git rev-parse @)
	REMOTE=$(git rev-parse "${UPSTREAM}")
	BASE=$(git merge-base @ "${UPSTREAM}")

	REPO_STATUS="\033[97;41;1m Diverged"
	if [ "${LOCAL}" = "${REMOTE}" ]; then
		REPO_STATUS="\033[32;1m Up-to-date"
	elif [ ${LOCAL} = ${BASE} ]; then
		REPO_STATUS="\033[33;1m Need to pull"
	elif [ ${REMOTE} = ${BASE} ]; then
		REPO_STATUS="\033[35;1m Need to push"
	fi

	printf "\n\tGit workspace status: %b \033[0m\n\n" "${REPO_STATUS}"
}

check_updates

snakemake --use-conda -p --cores "${CORES}"
