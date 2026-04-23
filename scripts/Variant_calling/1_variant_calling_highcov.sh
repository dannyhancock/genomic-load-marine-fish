#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N variant-calling
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=180gb
#PBS -l walltime=02:00:00:00

set -e

module purge
module load tools
module load perl/5.36.1
module load openjdk/22
module load bcftools/1.20
module load samtools/1.20
module load parallel/20220422

SPECIES="atlantic_cod"
THREADS_PER_CHROM=1
MAX_JOBS=23

BBC_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5_deduplicated_data/${SPECIES}"
VAR_OUT="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/6_snp_calling/${SPECIES}"
CHROM_OUT="$VAR_OUT/chromosomes_all"
mkdir -p "$VAR_OUT" "$CHROM_OUT"

REF="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}/genome.fna"
CHROMS="${VAR_OUT}/chromosomes.txt"

# Build BAM list
BAM_LIST_FILE="$VAR_OUT/all_bams.list"
: > "$BAM_LIST_FILE"

bams=( "$BBC_DIR"/*.rg.bam )

for bam in "${bams[@]}"; do
    [ -f "$bam" ] || continue
    base=$(basename "$bam")
    base=${base%.dedup_clipoverlap.bam}
    base=${base%.rg.bam}

    bai_default="${bam}.bai"
    bai_short="${bam%.bam}.bai"
    if [ ! -s "$bai_default" ] && [ ! -s "$bai_short" ]; then
        echo "Indexing $bam"
        samtools index -@ "$THREADS_PER_CHROM" "$bam"
    fi
    printf '%s\n' "$bam" >> "$BAM_LIST_FILE"
done

call_variants() {
    chrom=$1
    echo "  Processing chromosome: $chrom"
    bcftools mpileup \
        -r "$chrom" \
        -a FORMAT/DP,FORMAT/AD,FORMAT/SP \
        --fasta-ref "$REF" \
        --output-type u \
        --threads "$THREADS_PER_CHROM" \
        -b "$BAM_LIST_FILE" \
        --full-BAQ \
        -q 30 -Q 30 \
    | bcftools call \
        -a GQ,GP \
        --variants-only \
        --multiallelic-caller \
        --output-type u \
        --threads "$THREADS_PER_CHROM" \
    | bcftools view \
        --include 'QUAL>=30 && MQ>=30' \
        --output-type z \
        --output "$CHROM_OUT/${chrom}_variant.snps.bi.vcf.gz" \
        -m2 -M2 \
        --exclude-types indels \
        --threads "$THREADS_PER_CHROM"

    bcftools index -f --threads "$THREADS_PER_CHROM" "$CHROM_OUT/${chrom}_variant.snps.bi.vcf.gz"
}
export -f call_variants
export REF THREADS_PER_CHROM CHROM_OUT BAM_LIST_FILE

parallel --will-cite -j "$MAX_JOBS" call_variants ::: $(cat "$CHROMS")
wait
echo "Finished variant calling for all chromosomes"

# Merge per-chromosome VCFs into a single species-level VCF
MERGED="$VAR_OUT/${SPECIES}.QUAL30.MQ30.chroms.snps.bi.vcf.gz"
bcftools concat -a -Oz -o "$MERGED" "$CHROM_OUT/"*_variant.snps.bi.vcf.gz
bcftools index "$MERGED" --threads "$THREADS_PER_CHROM"

echo "Variant calling and concatenation complete: $MERGED"
