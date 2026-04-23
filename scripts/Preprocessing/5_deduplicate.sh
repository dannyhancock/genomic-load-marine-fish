#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N deduplicate-picard
#PBS -m n
#PBS -l nodes=1:ppn=20
#PBS -l mem=75gb
#PBS -l walltime=4:00:00

set -e

module purge
module load tools
module load openjdk/21
module load picard-tools/3.1.0

SPECIES="atlantic_cod"
DATA_TYPE="High_coverage_processed"

BASE="/home/projects/dp_00007/data/BlueBioClimate/Data/${DATA_TYPE}"
MAPPED_DIR="${BASE}/4_mapped_data/${SPECIES}"
OUT_DIR="${BASE}/5_deduplicated_data/${SPECIES}"
mkdir -p "$OUT_DIR"

for bam in "$MAPPED_DIR"/*.sorted.bam; do
    base=$(basename "$bam" .sorted.bam)

    echo "Deduplicating sample: $base"
    java -Xmx100G -jar "$PICARD" MarkDuplicates \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT \
        I="$bam" \
        O="$OUT_DIR/${base}.dedup.bam" \
        M="$OUT_DIR/${base}.duplication_metrics.txt" &
done
wait
echo "All samples deduplicated"


for bam in "$OUT_DIR"/*.dedup.bam; do
	base=$(basename "$bam" .dedup.bam)
	echo "Clipping overlap for sample: $base"
	/services/tools/bamutil/1.0.14/bam clipOverlap \
		--in "$bam" \
		--out "$OUT_DIR/${base}.dedup_clipoverlap.bam" \
		--stats > "$OUT_DIR/${base}.clipOverlap_stats.txt" &
done
wait
echo "All samples overlap clipped."

