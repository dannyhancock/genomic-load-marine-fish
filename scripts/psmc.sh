#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N psmc-atlantic_cod
#PBS -e /home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/21_demography/logs/psmc-atlantic_cod.err
#PBS -o /home/projects/dp_00007/data/BlueBioClimate/Data/High_coverage_processed/21_demography/logs/psmc-atlantic_cod.out
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=175gb
#PBS -l walltime=72:00:00

module load tools
module load psmc/0.6.5
module load bcftools/1.20
module load samtools/1.20
module load bedtools/2.31.1
module load gnuplot/6.0.0
module load qt/5.15.9

###############################################################################
# SETTINGS
###############################################################################
# for SPECIES in atlantic_cod atlantic_herring european_flounder lumpfish turbot; do
for SPECIES in atlantic_cod; do

	# # For plotting: scales axes from coalescent units to real time/Ne
	# # WE DO THIS IN R INSTEAD #
	# MU=1.64e-8      # mutation rate per site per generation
	# GEN_TIME=3     # generation time in years

	# Bootstrap and parallelization
	NR_BOOTSTRAPS=100
	MAX_JOBS=40

	# PSMC parameters
	# PSMC_PARAMS="-N25 -t15 -r5 -p 4+25*2+4+6"

	# Input paths
	BASE_DIR=/home/projects/dp_00007/data/BlueBioClimate/Data
	HC_DIR=${BASE_DIR}/High_coverage_processed
	BAM_DIR=${HC_DIR}/5_deduplicated_data/${SPECIES}

	# Reference genome
	REF=/home/projects/dp_00007/data/BlueBioClimate/Reference_genomes/friendly_name/${SPECIES}/genome.fna

	# Inversion regions to exclude: create complement BED for -R inclusion
	INVERSION_BED=${BASE_DIR}/00_info/${SPECIES}/inverted_regions.0based.bed
	MPILEUP_REGIONS=""
	HAS_INVERSIONS=false
	if [ -f "${INVERSION_BED}" ]; then
		echo "Will exclude inversions from: ${INVERSION_BED}"
		HAS_INVERSIONS=true
	fi

	# Output directory — use _noinv suffix when excluding inversions
	if [ "${HAS_INVERSIONS}" = true ]; then
		OUTDIR=${HC_DIR}/21_demography/${SPECIES}/psmc_noinv
	else
		OUTDIR=${HC_DIR}/21_demography/${SPECIES}/psmc
	fi
	mkdir -p ${OUTDIR}/consensus
	mkdir -p ${OUTDIR}/psmcfa
	mkdir -p ${OUTDIR}/results
	mkdir -p ${OUTDIR}/bootstrap
	mkdir -p ${OUTDIR}/plots

	# Generate keep-regions BED (complement of inversions) for bcftools -R
	if [ "${HAS_INVERSIONS}" = true ]; then
		KEEP_BED=${OUTDIR}/keep_regions.bed
		awk -v OFS='\t' '{print $1, $2}' ${REF}.fai > ${OUTDIR}/genome.txt
		bedtools complement -i ${INVERSION_BED} -g ${OUTDIR}/genome.txt > ${KEEP_BED}
		MPILEUP_REGIONS="-R ${KEEP_BED}"
	fi

	###############################################################################
	# STEP 1: Generate consensus FASTQ for each individual
	###############################################################################
	echo "Step 1: Generating consensus sequences..."

	# File to store mean depths per sample
	DEPTH_INFO_DIR=${BASE_DIR}/High_coverage_processed/00_info/${SPECIES}
	mkdir -p ${DEPTH_INFO_DIR}
	DEPTH_FILE=${DEPTH_INFO_DIR}/high_coverage_mean_depths.txt

	# Write header if file doesn't exist
	if [ ! -f ${DEPTH_FILE} ]; then
		echo -e "Sample\tMean_Depth" > ${DEPTH_FILE}
	fi

	for BAM in ${BAM_DIR}/*.rg.bam; do
		SAMPLE=$(basename ${BAM} .rg.bam)
		FQ_OUT=${OUTDIR}/consensus/${SAMPLE}.fq.gz

		if [ ! -f ${FQ_OUT} ]; then
			echo "Processing ${SAMPLE}..."

			# Calculate mean depth for this sample
			MEAN_DEPTH=$(samtools depth -a ${BAM} | awk '{sum+=$3; n++} END {print sum/n}')

			# Set MIN_DEPTH to 1/3 of mean, MAX_DEPTH to 2x mean (integer)
			MIN_DEPTH=$(awk "BEGIN {printf \"%d\", ${MEAN_DEPTH} / 3}")
			MAX_DEPTH=$(awk "BEGIN {printf \"%d\", ${MEAN_DEPTH} * 2}")

			echo "  Mean depth: ${MEAN_DEPTH}, using -d ${MIN_DEPTH} -D ${MAX_DEPTH}"

			# Record mean depth
			echo -e "${SAMPLE}\t${MEAN_DEPTH}" >> ${DEPTH_FILE}

			(
				bcftools mpileup -C50 ${MPILEUP_REGIONS} -f ${REF} ${BAM} | \
					bcftools call -c | \
					vcfutils.pl vcf2fq -d ${MIN_DEPTH} -D ${MAX_DEPTH} | \
					gzip > ${FQ_OUT}
			) &
		fi

		while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
			sleep 1
		done
	done

	wait
	echo "Consensus generation complete."

	###############################################################################
	# STEP 2: Convert FASTQ to PSMCFA format
	###############################################################################
	echo "Step 2: Converting to PSMCFA format..."

	for FQ in ${OUTDIR}/consensus/*.fq.gz; do
		SAMPLE=$(basename ${FQ} .fq.gz)
		PSMCFA_OUT=${OUTDIR}/psmcfa/${SAMPLE}.psmcfa

		if [ ! -f ${PSMCFA_OUT} ]; then
			echo "Converting ${SAMPLE}..."
			(fq2psmcfa -q20 ${FQ} > ${PSMCFA_OUT}) &
		fi

		while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
			sleep 1
		done
	done

	wait
	echo "PSMCFA conversion complete."

	###############################################################################
	# STEP 3: Run PSMC for each individual
	###############################################################################
	echo "Step 3: Running PSMC..."
	
	# Test combining more intervals
	TEST_PARAMS=(
		"4+25*2+4+6"
		"6+23*2+4+6"
		"2+2+25*2+4+6"
	)
	TEST_NAMES=(
		"4"
		"6"
		"2_2"
	)

	for PSMCFA in ${OUTDIR}/psmcfa/*.psmcfa; do
		SAMPLE=$(basename ${PSMCFA} .psmcfa)

		for i in "${!TEST_PARAMS[@]}"; do
			PSMC_PARAMS=${TEST_PARAMS[$i]}
			PARAM_ID=${TEST_NAMES[$i]}

			# Create output directory for this parameter set
			RESULTS_DIR=${OUTDIR}/results_${PARAM_ID}
			mkdir -p ${RESULTS_DIR}
			
			PSMC_OUT=${RESULTS_DIR}/${SAMPLE}.psmc

			echo "Running PSMC on ${SAMPLE} with params ${PSMC_PARAMS}..."
			psmc -N25 -t15 -r5 -p "${PSMC_PARAMS}" -o ${PSMC_OUT} ${PSMCFA} &

			while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
				sleep 1
			done
		done
	done

	wait
	echo "PSMC analysis complete."

	# ###############################################################################
	# # STEP 4: Bootstrap for confidence intervals
	# ###############################################################################
	# echo "Step 4: Running bootstrap replicates..."

	# for PSMCFA in ${OUTDIR}/psmcfa/*.psmcfa; do
		# SAMPLE=$(basename ${PSMCFA} .psmcfa)
		# BOOTSTRAP_DIR=${OUTDIR}/bootstrap/${SAMPLE}
		# mkdir -p ${BOOTSTRAP_DIR}

		# SPLIT_PSMCFA=${BOOTSTRAP_DIR}/${SAMPLE}.split.psmcfa
		# if [ ! -f ${SPLIT_PSMCFA} ]; then
			# echo "Splitting ${SAMPLE} for bootstrap..."
			# splitfa ${PSMCFA} > ${SPLIT_PSMCFA}
		# fi

		# echo "Running ${NR_BOOTSTRAPS} bootstrap replicates for ${SAMPLE}..."
		# seq ${NR_BOOTSTRAPS} | xargs -P ${MAX_JOBS} -I {} \
			# psmc ${PSMC_PARAMS} -b -o ${BOOTSTRAP_DIR}/${SAMPLE}.round-{}.psmc ${SPLIT_PSMCFA}
	# done

	# echo "Bootstrap complete."

	# # ###############################################################################
	# # # STEP 5: Generate plots
	# # ###############################################################################
	# echo "Step 5: Generating plots..."

	# for PSMC_RESULT in ${OUTDIR}/results/*.psmc; do
		# SAMPLE=$(basename ${PSMC_RESULT} .psmc)
		
		# echo "Plotting ${SAMPLE}..."
		# psmc_plot.pl -u ${MU} -g ${GEN_TIME} \
			# ${OUTDIR}/plots/${SAMPLE} \
			# ${PSMC_RESULT}
		
		# BOOTSTRAP_DIR=${OUTDIR}/bootstrap/${SAMPLE}

		# COMBINED=${OUTDIR}/plots/${SAMPLE}.combined.psmc
		# cat ${PSMC_RESULT} ${BOOTSTRAP_DIR}/${SAMPLE}.round-*.psmc > ${COMBINED}

		# echo "Plotting bootstrapped ${SAMPLE}..."
		# psmc_plot.pl -u ${MU} -g ${GEN_TIME} \
			# ${OUTDIR}/plots/${SAMPLE} \
			# ${COMBINED}
	# done

	# echo "Creating combined plot..."
	# psmc_plot.pl -u ${MU} -g ${GEN_TIME} \
		# ${OUTDIR}/plots/all_individuals \
		# ${OUTDIR}/results/*.psmc
		
	# # Convert to pdf
	# for f in ${OUTDIR}/plots/*.eps; do epstopdf "$f"; done

	# echo "Plots generated in ${OUTDIR}/plots/"

	###############################################################################
	# DONE
	###############################################################################
	echo "PSMC analysis complete for ${SPECIES}"
done
