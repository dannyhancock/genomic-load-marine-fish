#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N outgroup-variant-calling-atlantic_cod
#PBS -e /home/projects/dp_00007/data/BlueBioClimate/Data/Outgroups/6_variant_calling/logs/atlantic_cod_outgroup_variant_call.err
#PBS -o /home/projects/dp_00007/data/BlueBioClimate/Data/Outgroups/6_variant_calling/logs/atlantic_cod_outgroup_variant_call.out
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=150gb
#PBS -l walltime=3:00:00

set -e

module purge
module load tools
module load perl/5.36.1
module load openjdk/22
module load bcftools/1.20
module load samtools/1.20
module load parallel/20220422

# Variant calling and hard quality filtering

BAM_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Outgroups/5_deduplicated_data/atlantic_cod"
VAR_OUT="/home/projects/dp_00007/data/BlueBioClimate/Data/Outgroups/6_variant_calling/atlantic_cod"
mkdir -p "$VAR_OUT"

REF="/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/atlantic_cod/genome.fna"

# Variant sites found in focal species in which to find invariant sites
ORIG_REGIONS="/home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/6_snp_calling/atlantic_cod/atlantic_cod.QUAL30.MQ30.vcf.gz"
CHROMS="/home/projects/dp_00007/data/BlueBioClimate/Data/Outgroups/6_variant_calling/atlantic_cod/chromosomes.txt"

THREADS_PER_CHROM=4
TOTAL_THREADS=40
MAX_JOBS=$((TOTAL_THREADS / THREADS_PER_CHROM))

if [ ! -f "${ORIG_REGIONS}.csi" ]; then
	echo "Index for $ORIG_REGIONS not found, indexing..."
	bcftools index "$ORIG_REGIONS" --threads $THREADS_PER_CHROM
fi

# Pre-generate BED region files per chromosome from the known sites in the FOCAL species VCF.
for chrom in $(cat "$CHROMS"); do
    BED_FILE="$VAR_OUT/regions_${chrom}.bed"
    if [ ! -f "$BED_FILE" ]; then
        echo "Creating BED region file for $chrom..."
        bcftools query -f '%CHROM\t%POS0\t%POS\n' -r "$chrom" "$ORIG_REGIONS" > "$BED_FILE"
    fi
done

call_variants() {
    chrom=$1
    BED_FILE="$VAR_OUT/regions_${chrom}.bed"
    echo "  Processing chromosome: $chrom using region file: $BED_FILE"
    
	# Find all variant sites on this chromosome
    bcftools mpileup \
		-r $chrom \
		-a FORMAT/DP \
		--fasta-ref "$REF" \
        --output-type u \
        --threads $THREADS_PER_CHROM \
        "$bam" \
    | bcftools call \
		-a GQ,GP \
		--variants-only \
        --multiallelic-caller \
        --output-type u \
        --threads $THREADS_PER_CHROM \
    | bcftools view \
        --include 'QUAL>=30 && MQ>=30' \
        --output-type z \
        --output "$SAMPLE_OUT/${sample}_${chrom}_variant.vcf.gz" \
        --threads $THREADS_PER_CHROM
		
	bcftools index --threads $THREADS_PER_CHROM "$SAMPLE_OUT/${sample}_${chrom}_variant.vcf.gz"
	
	# Repeat to find invariant sites at sites where variants were found in Focal species
	bcftools mpileup \
        -r $chrom \
		-T "$BED_FILE" \
		-a FORMAT/DP \
		--fasta-ref "$REF" \
        --output-type u \
        --threads $THREADS_PER_CHROM \
        "$bam" \
    | bcftools call \
		-a GQ,GP \
        --multiallelic-caller \
        --output-type u \
        --threads $THREADS_PER_CHROM \
    | bcftools view \
        --include 'QUAL>=30 && MQ>=30 && ALT=="." && INFO/INDEL!=1' \
        --output-type z \
        --output "$SAMPLE_OUT/${sample}_${chrom}_invariant.vcf.gz" \
        --threads $THREADS_PER_CHROM
		
	bcftools index --threads $THREADS_PER_CHROM "$SAMPLE_OUT/${sample}_${chrom}_invariant.vcf.gz"
}
export -f call_variants

# Call variants on chromosomes in parallel for one outgroup sample after another
for bam in "$BAM_DIR"/*.dedup_clipoverlap.bam; do
    sample=$(basename "$bam" .dedup_clipoverlap.bam)
	
    echo "Variant calling for sample: $sample"
    
    # Check for BAM index; if not present, create it.
    if [ ! -f "${bam}.bai" ]; then
         echo "Index for $bam not found, indexing..."
         samtools index "$bam"
    fi

    # Create an output directory for per-chromosome VCFs for this sample.
    SAMPLE_OUT="$VAR_OUT/${sample}"
    mkdir -p "$SAMPLE_OUT"
	
	# Export variables required by call_variants.
	export bam sample SAMPLE_OUT VAR_OUT REF THREADS_PER_CHROM
	
	# Run variant calling for each chromosome in parallel
	parallel --will-cite -j $MAX_JOBS call_variants ::: $(cat "$CHROMS")
    
    # Wait for all chromosome jobs for this sample to complete.
    wait
    echo "  Finished variant calling for all chromosomes of sample: $sample"
    
    # Merge per-chromosome VCFs into a single VCF file for the sample.
    bcftools concat -a -O z -o "$VAR_OUT/${sample}_variant.vcf.gz" "$SAMPLE_OUT/"*_variant.vcf.gz
	bcftools index "$VAR_OUT/${sample}_variant.vcf.gz" --threads $THREADS_PER_CHROM
	# Merge per-chromosome INVARIANT sites into a single VCF for the sample
	bcftools concat -a -O z -o "$VAR_OUT/${sample}_invariant.vcf.gz" "$SAMPLE_OUT/"*_invariant.vcf.gz
	bcftools index "$VAR_OUT/${sample}_invariant.vcf.gz" --threads $THREADS_PER_CHROM
    
	# Create a single VCF for the sample and sort (using a temp file) and index
	bcftools concat -a -O z -o "$VAR_OUT/${sample}.vcf.gz" "$VAR_OUT/${sample}_variant.vcf.gz" "$VAR_OUT/${sample}_invariant.vcf.gz"
	bcftools sort -O z -o "$VAR_OUT/${sample}.vcf.sorted.gz" "$VAR_OUT/${sample}.vcf.gz"
	mv "$VAR_OUT/${sample}.vcf.sorted.gz" "$VAR_OUT/${sample}.vcf.gz"
	bcftools index "$VAR_OUT/${sample}.vcf.gz" --threads $THREADS_PER_CHROM
	
	# remove indels and multiallelic sites (keep invariant sites with m1)
	bcftools view "$VAR_OUT/${sample}.vcf.gz" --threads $THREADS_PER_CHROM \
	-m1 -M2 --exclude-types indels -Oz -o "$VAR_OUT/${sample}.noindels.bi.vcf.gz"
	bcftools index "$VAR_OUT/${sample}.noindels.bi.vcf.gz" --threads $THREADS_PER_CHROM
	
	echo "Variant calling and concatenation complete for sample: $sample"
done

echo "Variant calling for all samples completed."