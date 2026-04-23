#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N fastqc-trimmed
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=175gb
#PBS -l walltime=1:00:00

set -e

module purge
module load tools
module load perl/5.36.1
module load openjdk/21
module load fastqc/0.11.9
module load multiqc/1.12

SPECIES="atlantic_cod"
DATA_TYPE="High_coverage_processed"

BASE="/home/projects/dp_00007/data/BlueBioClimate/Data/${DATA_TYPE}"
TRIMDIR="${BASE}/2_trimmed_data/${SPECIES}"
FAST_OUT="${BASE}/3_qc_trimmed/${SPECIES}/fastqc"
MULTI_OUT="${BASE}/3_qc_trimmed/${SPECIES}/multiqc"

mkdir --parents "$FAST_OUT" "$MULTI_OUT"

trimmed_files=$(find "$TRIMDIR" -type f -name '*fq.gz')

fastqc -t 40 -q -o "${FAST_OUT}" ${trimmed_files}
echo "FastQC reports generated"

multiqc --force \
    --filename "${MULTI_OUT}/multiqc_${SPECIES}_${DATA_TYPE}_trimmed.html" \
    "$FAST_OUT"

echo "MultiQC report generated"
