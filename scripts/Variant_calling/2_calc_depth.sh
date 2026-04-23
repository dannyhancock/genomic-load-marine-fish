#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N depth_all_species
#PBS -e /home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/0_depth/logs/depth_all_species.err
#PBS -o /home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/0_depth/logs/depth_all_species.out
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=180gb
#PBS -l walltime=48:00:00

set -e

module purge
module load tools
module load htslib/1.20
module load angsd/0.940

THREADS_PER_SPECIES=8

# Samples to drop (QC failures, library controls, duplicates) from all 5 species
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

run_species() {
    local SPECIES="$1"
    local GENOME

    case "$SPECIES" in
        atlantic_cod)
            GENOME="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}" ;;
        atlantic_herring)
            GENOME="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}" ;;
        european_flounder)
            GENOME="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}" ;;
        lumpfish)
            GENOME="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}" ;;
        turbot)
            GENOME="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}" ;;
        *)
            echo "ERROR: unknown species $SPECIES" >&2
            return 5 ;;
    esac

    local OUTROOT="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/0_depth/${SPECIES}"
    mkdir -p "$OUTROOT"

    # Build BAM filelist: low-coverage + downsampled high-coverage
    local filelist="${OUTROOT}/bam.filelist"
    : > "$filelist"
    realpath --quiet "/home/projects/dp_00007/data/BlueBioClimate/Data/Low_coverage_processed/5_deduplicated_data/${SPECIES}"/*.rg.bam \
        >> "$filelist"
    realpath --quiet "/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5.2_downsampling/${SPECIES}"/*.rg.bam \
        >> "$filelist"

    for sample in "${exclude_samples[@]}"; do
        sed -i "/$sample/d" "$filelist"
    done

    local NIND
    NIND=$(wc -l < "$filelist")
    local MIN_IND=$(( NIND / 2 ))
    echo "[$SPECIES] NIND=$NIND  MIN_IND=$MIN_IND" >&2

    angsd -bam "$filelist" \
        -ref "$GENOME/genome.fna" \
        -fai "$GENOME/genome.fna.fai" \
        -nThreads "$THREADS_PER_SPECIES" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -trim 0 \
        -C 50 \
        -baq 2 \
        -minQ 20 \
        -minMapQ 20 \
        -doDepth 1 \
        -doCounts 1 \
        -maxDepth 2000 \
        -minInd "$MIN_IND" \
        -out "${OUTROOT}/${SPECIES}_lc_and_downsampled_realigned"

    echo "[$SPECIES] Depth distribution written to: ${OUTROOT}/${SPECIES}_lc_and_downsampled_realigned.*"
}
export -f run_species
export THREADS_PER_SPECIES exclude_samples

SPECIES_LIST=(atlantic_cod atlantic_herring european_flounder lumpfish turbot)

for SPECIES in "${SPECIES_LIST[@]}"; do
    run_species "$SPECIES" &
done
wait

echo "Depth distributions complete for all species."
