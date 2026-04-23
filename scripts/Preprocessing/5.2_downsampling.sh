#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N downsampling_atlantic_cod
#PBS -e /home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5.2_downsampling/logs/downsampling_atlantic_cod.err
#PBS -o /home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5.2_downsampling/logs/downsampling_atlantic_cod.out
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=150gb
#PBS -l walltime=00:30:00

set -e

module purge
module load tools
module load perl/5.36.1
module load vcftools/0.1.16
module load samtools/1.20
module load parallel/20220422

## DIRECTORIES ##
HIGHCOV_SNPS_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/6_snp_calling/atlantic_cod"
LOWCOV_SNPS_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Low_coverage_processed/6_snp_calling/atlantic_cod"

HIGHCOV_BAM_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5_deduplicated_data/atlantic_cod"
LOWCOV_BAM_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Low_coverage_processed/5_deduplicated_data/atlantic_cod"

DOWNSAMPLE_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/5.2_downsampling/atlantic_cod"
mkdir -p $DOWNSAMPLE_DIR

#################################################
## Calculate depth of the low coverage samples ##
#################################################
# Output file for mean depths of low-coverage samples
lowcov_mean_file="$DOWNSAMPLE_DIR/low_coverage_mean_depths.txt" 
> "$lowcov_mean_file"

# Calculate the mean depth for each low-coverage BAM file
for bamfile in "$LOWCOV_BAM_DIR"/*.rg.bam; do
    sample_name=$(basename "$bamfile" .rg.bam)

    # Calculate mean depth for each BAM file, including sites with 0 depth
    mean_depth=$(samtools depth -a "$bamfile" | awk '{sum += $3} END {if (NR > 0) print sum / NR}')

    # Append sample mean depth to the output file
    echo "$sample_name: $mean_depth" >> "$lowcov_mean_file"
done

# Calculate the overall mean depth for all low-coverage BAMs
target_depth=$(awk '{sum += $2; count++} END {if (count > 0) print sum / count}' "$lowcov_mean_file")

echo $target_depth > $DOWNSAMPLE_DIR/low_coverage_mean.txt

echo "Target depth for downsampling high coverage BAMs: $target_depth"

##################
## DOWNSAMPLING ##
##################
# Output file for mean depths of all high coverage samples
highcov_mean_file="$DOWNSAMPLE_DIR/high_coverage_mean_depths.txt" 
> "$highcov_mean_file"

downsampled_mean_file="$DOWNSAMPLE_DIR/downsampled_seeds_depths.txt"
echo -e "Sample\tMean_depth\tSeed" > "$downsampled_mean_file"

# Generate a single file list with all high-coverage BAM files
filelist="$DOWNSAMPLE_DIR/atlantic_cod_highcov.filelist"
ls "$HIGHCOV_BAM_DIR"/*.rg.bam > "$filelist"

# Remove high coverage samples that were also sampled at low depth
exclude_samples=(
    "OOB07_135_EKDL230051351-1A_HVLNVDSX7_L2"
    "OOB07_136_EKDL230051351-1A_HVLNVDSX7_L2"
    "OOB07_142_EKDL230051351-1A_HVLNVDSX7_L2"
    "OOB07_149_EKDL230051351-1A_HVLNVDSX7_L2"
    "OOB07_151_EKDL230051351-1A_HVLNVDSX7_L2"
    "OOB07_152_EKDL230051351-1A_HVLNVDSX7_L2"
    "Kat15_02_EKDL230051351-1A_HVLNVDSX7_L2"
    "Kat15_08_EKDL230051351-1A_HVLNVDSX7_L2"
    "Kat15_10_EKDL230051351-1A_HVLNVDSX7_L2"
    "Kat15_13_EKDL230051351-1A_HVLNVDSX7_L2"
    "Kat15_14_EKDL230051351-1A_HVLNVDSX7_L2"
    "Kat15_17_EKDL230051351-1A_HVLNVDSX7_L2"
    "MorayFirth2003_09_EKDL230051351-1A_HVLNVDSX7_L2"
    "MorayFirth2003_16_EKDL230051351-1A_HVLNVDSX7_L2"
    "MorayFirth2003_22_EKDL230051351-1A_HVLNVDSX7_L2"
    "MorayFirth2003_35_EKDL230051351-1A_HVLNVDSX7_L2"
    "MorayFirth2003_38_EKDL230051351-1A_HVLNVDSX7_L2"
    "Tur_29_9011_EKDL220000316-1a-AK31123-AK31124_HWHMNDSX2_L4"
)
for sample in "${exclude_samples[@]}"; do
    sed -i "/$sample/d" "$filelist"
done

process_bam() {
    bamfile=$1
    sample_name=$(basename "$bamfile" .rg.bam)
    tmp_bam="$DOWNSAMPLE_DIR/${sample_name}_tmp.bam"  # Temporary BAM during processing
    output_bam="$DOWNSAMPLE_DIR/${sample_name}_downsampled.rg.bam"  # Final BAM after processing

    # Calculate mean depth for the current high-coverage BAM
    mean_depth=$(samtools depth -a "$bamfile" | awk '{sum += $3} END {if (NR > 0) print sum / NR}')
    if (( $(echo "$mean_depth == 0" | bc -l) )); then
        echo "Warning: $sample_name has zero depth, skipping..."
        return
    fi

    # Calculate the downsampling fraction
    subsample_fraction=$(echo "$target_depth / $mean_depth" | bc -l)

    # Generate a numeric seed from the sample name
    seed=$(echo -n "$sample_name" | od -An -t u4 | awk '{s += $1} END {print s % 10000}')

    # Downsample BAM to a temporary file
    samtools view -bh --subsample-seed "$seed" --subsample "$subsample_fraction" "$bamfile" -o "$tmp_bam" --threads $THREADS_PER_BAM
    if [[ ! -f "$tmp_bam" ]]; then
        echo "Error: Temporary downsampled BAM for $sample_name not created, skipping..."
        return
    fi

    # Add or replace read group to create the final output BAM
    samtools addreplacerg -r "ID:${sample_name}_downsampled" -r "SM:${sample_name}_downsampled" -r "PL:Illumina" -o "$output_bam" "$tmp_bam"
     # Clean up the intermediate BAM file
	rm "$tmp_bam"
	
	samtools index "$output_bam"

    # Calculate depth stats for the final BAM
    downsampled_mean_depth=$(samtools depth -a "$output_bam" | awk '{sum += $3} END {if (NR > 0) print sum / NR}')
    echo -e "$sample_name\t$downsampled_mean_depth\t$seed" >> "$downsampled_mean_file"

    echo "Downsampled $sample_name to approximately target depth $target_depth"
}

export -f process_bam

THREADS_PER_BAM=4 
TOTAL_THREADS=40
MAX_JOBS=$((TOTAL_THREADS / THREADS_PER_BAM))

# Export variables for parallel
export DOWNSAMPLE_DIR target_depth highcov_mean_file downsampled_mean_file THREADS_PER_BAM

# Run in parallel
cat "$filelist" | parallel --will-cite -j $MAX_JOBS process_bam {}


