#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N fastp-trim
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=150gb
#PBS -l walltime=4:00:00

set -e

module purge
module load tools
module load fastp/0.23.2

SPECIES="atlantic_cod"
DATA_TYPE="High_coverage_processed"
THREADS_PER_JOB=10

BASE="/home/projects/dp_00007/data/BlueBioClimate/Data/${DATA_TYPE}"
RAWDIR="${BASE}/0_raw_data/${SPECIES}"
TRIM_OUT="${BASE}/2_trimmed_data/${SPECIES}"
mkdir -p "$TRIM_OUT"

tmpfile=$(mktemp)
find "$RAWDIR" -type f \( -name '*_1.fq.gz' -o -name '*_1.fastq.gz' \) > "$tmpfile"

while read -r in1; do
    in2=$(echo "$in1" | sed 's/_1\./_2./')
    if [ ! -f "$in2" ]; then
        echo "Warning: paired file for $in1 not found ($in2). Skipping."
        continue
    fi

    base=$(basename "$in1")
    base=${base/_1.fq.gz/}
    base=${base/_1.fastq.gz/}

    out1="$TRIM_OUT/${base}_1.fq.gz"
    out2="$TRIM_OUT/${base}_2.fq.gz"

    echo "Trimming sample: $base"
    fastp \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --trim_poly_g \
        --thread "$THREADS_PER_JOB" \
        --length_required 35 \
        --low_complexity_filter \
        --complexity_threshold 30 \
        --correction \
        --in1 "$in1" \
        --in2 "$in2" \
        --out1 "$out1" \
        --out2 "$out2" &
done < "$tmpfile"

wait
echo "All samples trimmed"
rm "$tmpfile"
