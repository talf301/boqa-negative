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
data=/dupa-filer/talf/boqa-negative/

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

logdir=~/sge_logs/boqa/`basename $out`
mkdir -pv "$logdir"
mkdir -pv "$out"
mkdir -pv "$out/scripts"

for file in $dir/*_hpo.txt; do
	f=`basename $file _hpo.txt`
	if [[ -s "$out/$f"_hpo.txt.rank ]]; then
		echo "Output file already exists "$out/$f".finished" >&2
		continue
	fi
	script="$out/scripts/dispatch_$f.sh"
	cat > "$script" <<EOF
#!/usr/bin/env bash
#$ -V
#$ -N "BOQA_$f"
#$ -l h_vmem="$memory"
#$ -e $logdir
#$ -o $logdir
#$ -l hostname="supa*"

set -eu
set -o pipefail
temp=\$TMPDIR/$f

mkdir -pv \$temp

echo "Current directory: \$(pwd)" >&2
echo "Temp directory: \$temp" >&2
echo "Input file: $dir/$f"_hpo.txt >&2
echo "Target output director: $out" >&2

test -s "$dir/$f"_hpo.txt

# Run the actual python script
/filer/tools/python/Python-2.7.6/python $data/run_net.py -D $data -P "$dir/$f"_hpo.txt -O \$temp -k 10

# Make sure moving finishes correctly
mv -v \$temp/* $out
touch "$out/$f".finished
EOF

	# Wait for space on cluster
	sge_wait
	# Submit job
	qsub -S /bin/sh "$script"
done
