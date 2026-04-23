#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N bwa-mem
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=150gb
#PBS -l walltime=24:00:00

set -e

module purge
module load tools
module load bwa/0.7.17
module load samtools/1.18
module load gnuplot/6.0.0

SPECIES="atlantic_cod"
DATA_TYPE="High_coverage_processed"
NCPU=10
MAX_JOBS=4

BASE="/home/projects/dp_00007/data/BlueBioClimate/Data/${DATA_TYPE}"
TRIMDIR="${BASE}/2_trimmed_data/${SPECIES}"
MAP_OUT="${BASE}/4_mapped_data/${SPECIES}"
INDEX="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/4_mapped_data/bwa_indexes/${SPECIES}"
SCRATCH="/home/projects/dp_00007/scratch"

mkdir -p "$MAP_OUT"

for in1 in "$TRIMDIR"/*_1.fq.gz; do
    base=$(basename "$in1")
    base=${base/_1.fq.gz/}
    in2="$TRIMDIR/${base}_2.fq.gz"

    if [ ! -f "$in2" ]; then
        echo "Paired file for $in1 not found. Skipping $base."
        continue
    fi

    {
        echo "Mapping sample: $base"
        RG="@RG\tID:$base\tSM:$base\tPL:Illumina"

        bwa mem -t "$NCPU" -M -R "$RG" "$INDEX" "$in1" "$in2" > "$MAP_OUT/${base}.sam"

        samtools stats "$MAP_OUT/${base}.sam" > "$MAP_OUT/${base}.stats"
        plot-bamstats "$MAP_OUT/${base}.stats" -p "$MAP_OUT/${base}_plots/"

        samtools view -bS -h -q 20 -F 4 -@ "$NCPU" "$MAP_OUT/${base}.sam" > "$MAP_OUT/${base}.bam"
        samtools sort -@ "$NCPU" -m 3G -T "$SCRATCH" -o "$MAP_OUT/${base}.sorted.bam" "$MAP_OUT/${base}.bam"
        samtools index -@ "$NCPU" "$MAP_OUT/${base}.sorted.bam"

        samtools stats "$MAP_OUT/${base}.sorted.bam" > "$MAP_OUT/${base}.filtered.stats"
        plot-bamstats "$MAP_OUT/${base}.filtered.stats" -p "$MAP_OUT/${base}_filtered_plots/"

        rm "$MAP_OUT/${base}.sam" "$MAP_OUT/${base}.bam"
        echo "Mapping completed for sample: $base"
    } &

    while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 5
    done
done

wait
echo "All samples processed."
