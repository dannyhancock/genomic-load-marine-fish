#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N ANGSD-SNP-downsampled
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=180gb
#PBS -l walltime=72:00:00

set -e

module purge
module load tools
module load bcftools/1.20
module load htslib/1.20
module load angsd/0.940
module load pcangsd/20220330
module load parallel/20220422

SPECIES="atlantic_cod"
REF="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/gadus_morhua_GCF_902167405.1/genome.fna"
ANGSD_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/1_ANGSD/${SPECIES}"
CHROM_DIR="${ANGSD_DIR}/chromosomes"
mkdir -p "$ANGSD_DIR" "$CHROM_DIR"

CHROMOSOMES="${ANGSD_DIR}/chromosomes.txt"

# Build BAM filelist (low coverage + downsampled high coverage)
filelist="/home/projects/dp_00007/data/BlueBioClimate/Data/Low_coverage_processed/6_snp_calling/file_lists/${SPECIES}/lc_and_downsampled.filelist"
mkdir -p "$(dirname "$filelist")"
ls "/home/projects/dp_00007/data/BlueBioClimate/Data/Low_coverage_processed/5_deduplicated_data/${SPECIES}"/*.rg.bam > "$filelist"
ls "/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5.2_downsampling/${SPECIES}"/*.rg.bam >> "$filelist"

# Drop QC failures, library controls and duplicate samples
exclude_samples=(
    "BBassin_2007_145_EKDL210007214-1a-AK11349-AK9138_HK5MFDSX2_L1"
    "Thy_04_22_EKDL210007215-1a-AK11420-AK18623_HK5MFDSX2_L2"
    "S19_18_EKDL210007215-1a-5UDI1665-5UDI1473_HK5MFDSX2_L2"
    "Tur_21_1739_EKDL210007215-1a-AK17448-AK31259_HK5MFDSX2_L2"
    "BB_2007_130A_02_EKDL210007214-1a-AK12650-AK23830_HK5MFDSX2_L1"
    "BB_2007_130A_1_EKDL210007213-1a-AK11653-AK16291_HK53LDSX2_L3"
    "BB_2007_130_02_EKDL210007214-1a-AK11322-AK30698_HK5MFDSX2_L1"
    "BB_2007_130_1_EKDL210007213-1a-AK11473-AK16281_HK53LDSX2_L3"
    "BBassin_2007_130A_EKDL210007215-1a-AK17466-AK22037_HK5MFDSX2_L2"
    "Flo_21_1113_01_EKDL210007213-1a-AK11387-AK16228_HK53LDSX2_L3"
    "Flo_21_1113_02_EKDL210007214-1a-AK17103-AK17021_HK5MFDSX2_L1"
    "Flo_21_1113A_01_EKDL210007213-1a-AK11615-AK16297_HK53LDSX2_L3"
    "Flo_21_1113A_02_EKDL210007214-1a-AK2228-AK30697_HK5MFDSX2_L1"
    "Flo_21_1113A_EKDL210007215-1a-AK17469-AK33136_HK5MFDSX2_L2"
    "Flo_21_1113_EKDL210007215-1a-AK11627-AK31252_HK5MFDSX2_L2"
    "Tur_29_9009_EKDL210007213-1a-AK11467-AK14066_HK53LDSX2_L3"
)
for sample in "${exclude_samples[@]}"; do
    sed -i "/$sample/d" "$filelist"
done

# Per-species depth limits
depth_limits="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/0_depth/depth_limits.txt"
minDepth=$(awk -v s="$SPECIES" '$1 == s {print $2}' "$depth_limits")
maxDepth=$(awk -v s="$SPECIES" '$1 == s {print $3}' "$depth_limits")

n_individuals=$(wc -l < "$filelist")
MIN_IND=$(( n_individuals / 2 ))

PREFIX="${SPECIES}_lc_and_downsampled.minDP${minDepth}.maxDP${maxDepth}.minInd${MIN_IND}"

THREADS_PER_CHROM=8
TOTAL_THREADS=40
MAX_JOBS=$(( TOTAL_THREADS / THREADS_PER_CHROM ))

# Per-chromosome ANGSD
process_chromosome() {
    chrom=$1
    out="${CHROM_DIR}/${PREFIX}.${chrom}"
    echo "Processing chromosome: ${chrom}"

    angsd -bam "$filelist" \
        -ref "${REF}" \
        -fai "${REF}.fai" \
        -nThreads "${THREADS_PER_CHROM}" \
        -r "${chrom}" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -trim 0 \
        -C 50 \
        -baq 1 \
        -minQ 20 \
        -minMapQ 20 \
        -setMinDepth "${minDepth}" \
        -setMaxDepth "${maxDepth}" \
        -nInd "${n_individuals}" \
        -GL 1 \
        -doMajorMinor 1 \
        -doGlf 2 \
        -doGeno 1 \
        -doPost 1 \
        -doMaf 1 \
        -SNP_pval "1e-6" \
        -doCounts 1 \
        -doDepth 1 \
        -doBcf 1 \
        -minInd "${MIN_IND}" \
        -out "${out}"

    echo "Chromosome ${chrom} done."
}
export -f process_chromosome
export REF CHROM_DIR PREFIX filelist minDepth maxDepth n_individuals MIN_IND THREADS_PER_CHROM

cat "${CHROMOSOMES}" | parallel --will-cite -j "${MAX_JOBS}" process_chromosome {}

#########################################################
# Stitch per-chromosome outputs into genome-wide files
#########################################################
echo "Concatenating per-chromosome outputs..."

# Order chromosomes as listed in CHROMOSOMES file
mapfile -t CHROM_ORDER < "${CHROMOSOMES}"

# Beagle: keep header from the first chrom, drop headers from the rest
first=1
for chrom in "${CHROM_ORDER[@]}"; do
    f="${CHROM_DIR}/${PREFIX}.${chrom}.beagle.gz"
    if [ "$first" -eq 1 ]; then
        zcat "$f"
        first=0
    else
        zcat "$f" | tail -n +2
    fi
done | gzip > "${ANGSD_DIR}/${PREFIX}.beagle.gz"

# MAFs: same header strategy
first=1
for chrom in "${CHROM_ORDER[@]}"; do
    f="${CHROM_DIR}/${PREFIX}.${chrom}.mafs.gz"
    if [ "$first" -eq 1 ]; then
        zcat "$f"
        first=0
    else
        zcat "$f" | tail -n +2
    fi
done | gzip > "${ANGSD_DIR}/${PREFIX}.mafs.gz"

# BCF: bcftools concat handles headers correctly
bcf_list="${CHROM_DIR}/${PREFIX}.bcf.list"
> "$bcf_list"
for chrom in "${CHROM_ORDER[@]}"; do
    bcf="${CHROM_DIR}/${PREFIX}.${chrom}.bcf"
    bcftools index -f "$bcf"
    echo "$bcf" >> "$bcf_list"
done
bcftools concat -f "$bcf_list" -Ob -o "${ANGSD_DIR}/${PREFIX}.bcf" --threads 8
bcftools index "${ANGSD_DIR}/${PREFIX}.bcf"

echo "Stitched genome-wide outputs:"
echo "  ${ANGSD_DIR}/${PREFIX}.beagle.gz"
echo "  ${ANGSD_DIR}/${PREFIX}.mafs.gz"
echo "  ${ANGSD_DIR}/${PREFIX}.bcf"

#####################################
## PCAngsd on the genome-wide beagle
#####################################
PCA_OUT="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/3_PCANGSD/${SPECIES}"
mkdir -p "$PCA_OUT"
pcangsd \
    --beagle "${ANGSD_DIR}/${PREFIX}.beagle.gz" \
    --out "${PCA_OUT}/${PREFIX}.PCAngsd" \
    --selection \
    --selection_e 2 \
    --threads 8

echo "GLs and SNPs called for low coverage and downsampled individuals"
