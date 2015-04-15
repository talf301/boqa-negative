#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
	cat <<EOF
Usage: $0 dir out

Dispatch boqa with sampling technique on all of the _hpo.txt files
in a directory.
EOF
	exit 1
}

if [ $# -ne 2 ]; then
	usage
fi

dir=$1
out=$2
max_jobs=100
memory=5g
processors=1

function sge_wait {
    sleep_time=1  # seconds
    # Check SGE for number of jobs
    local sge_jobs="$(qstat | grep "BOQA_" | wc -l)" || true
    while [[ $sge_jobs -ge $max_jobs ]]; do
	sleep $sleep_time
	sleep_time=$(expr $sleep_time "*" 2)
	sge_jobs="$(qstat | grep "BOQA_" | wc -l)" || true
    done
}

logdir=~/sge_logs/boqa/"$out"
mkdir -pv "$logdir"
mkdir -pv "$out"
mkdir -pv "$out/scripts"

for file in $dir/*_hpo.txt; do
	f=`basename $file _hpo.txt`
	
