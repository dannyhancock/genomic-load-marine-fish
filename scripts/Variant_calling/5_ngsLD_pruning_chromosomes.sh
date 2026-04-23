#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N ngsLD-pruning-chromosomes-atlantic_cod
#PBS -e /home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/2_ngsLD/logs/ngsLD_pruning_chromosomes_atlantic_cod.err
#PBS -o /home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/2_ngsLD/logs/ngsLD_pruning_chromosomes_atlantic_cod.out
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=180gb
#PBS -l walltime=40:00:00

module purge
module load tools
module load parallel/20220422
module load gsl/2.7

ANGSD_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/1_ANGSD/atlantic_cod"
ngsLD_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/2_ngsLD/atlantic_cod"
mkdir --parents $ngsLD_DIR
mkdir --parents $ngsLD_DIR/stripped # folder for storing stripped beagle/pos files for ngsLD

# # Get min and max DP from distributions
depth_limits="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/0_depth/depth_limits.txt"
species="atlantic_cod"
minDepth=$(awk -v species="$species" '$1 == species {print $2}' "$depth_limits")
maxDepth=$(awk -v species="$species" '$1 == species {print $3}' "$depth_limits")

filelist="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/00_cluster_lists/atlantic_cod/lc_and_downsampled.filelist"
n_individuals=$(wc --lines "$filelist" | cut -d " " -f 1)
MIN_IND=$(( n_individuals / 2 ))

# Filtered BEAGLE prefix (Filtered on min and max depth)
FILTERED_PREFIX="atlantic_cod_lc_and_downsampled.minDP${minDepth}.maxDP${maxDepth}.minInd${MIN_IND}"

# Software
ngsLD="/home/projects/dp_00007/data/BlueBioClimate/Data/Software/ngsLD"
prune_graph="/home/projects/dp_00007/data/BlueBioClimate/Data/Software/prune_graph"

CHROMOSOMES_FILE="$ANGSD_DIR/chromosomes.txt"
THREADS_PER_CHROM=4
TOTAL_THREADS=40
MAX_JOBS=$((TOTAL_THREADS / THREADS_PER_CHROM))

MAX_DIST=10000
MAX_KB_DIST=$(($MAX_DIST/1000))
R2=0.1

# Function to process each chromosome
process_chromosome() {
  chrom=$1
  
  # Skip if already done (non-empty file)
  if [[ -s "$ngsLD_DIR/${FILTERED_PREFIX}_${chrom}.ld" ]]; then
    echo "Skipping finding LD for $chrom - found $ngsLD_DIR/${FILTERED_PREFIX}_${chrom}.ld"
  else
  	  echo "Processing chromosome: $chrom"

	  # Prepare the Beagle file for the chromosome
	  zcat $ANGSD_DIR/$FILTERED_PREFIX.beagle.gz | awk -v chrom="$chrom" 'NR==1 || $1 ~ "^"chrom"_"' | tail -n +2 | cut -f 4- | gzip > $ngsLD_DIR/stripped/${FILTERED_PREFIX}_${chrom}_stripped.beagle.gz

	  # Prepare the positions file for the chromosome
	  zcat $ANGSD_DIR/$FILTERED_PREFIX.mafs.gz | awk -v chrom="$chrom" 'NR==1 || $1 == chrom' | tail -n +2 | cut -f 1,2 | sed 's/:/_/g' | gzip > $ngsLD_DIR/stripped/${FILTERED_PREFIX}_${chrom}_stripped.pos.gz

	  # Get the number of individuals and sites
	  n_inds=$(zcat $ngsLD_DIR/stripped/${FILTERED_PREFIX}_${chrom}_stripped.beagle.gz | head -n 1 | awk '{print NF/3}')
	  n_sites=$(zcat $ngsLD_DIR/stripped/${FILTERED_PREFIX}_${chrom}_stripped.beagle.gz | wc -l)

	  # Run ngsLD with N threads per chromosome
	  $ngsLD/ngsLD \
	  --geno $ngsLD_DIR/stripped/${FILTERED_PREFIX}_${chrom}_stripped.beagle.gz \
	  --pos $ngsLD_DIR/stripped/${FILTERED_PREFIX}_${chrom}_stripped.pos.gz \
	  --probs \
	  --n_ind $n_inds \
	  --n_sites $n_sites \
	  --max_kb_dist $MAX_KB_DIST \
	  --min_maf 0.05 \
	  --n_threads $THREADS_PER_CHROM \
	  --out "$ngsLD_DIR/${FILTERED_PREFIX}_${chrom}.ld"
  fi
	
  # Skip if already done (non-empty file)
  if [[ -s "$ngsLD_DIR/${FILTERED_PREFIX}_${chrom}_unlinked.pos" ]]; then
    echo "Skipping pruning for $chrom - found $ngsLD_DIR/${FILTERED_PREFIX}_${chrom}_unlinked.pos"
  else
	  # Prune linked SNPs
	  $prune_graph/target/release/prune_graph --header --in "$ngsLD_DIR/${FILTERED_PREFIX}_${chrom}.ld" \
	  --weight-field "r2" \
	  --weight-filter "dist <= $MAX_DIST && r2 >= ${R2}" \
	  -n $THREADS_PER_CHROM \
	  --out "$ngsLD_DIR/${FILTERED_PREFIX}_${chrom}_unlinked.pos" \
	  --verbose
  fi
  
  echo "Chromosome $chrom processed."
}

export -f process_chromosome

export ANGSD_DIR ngsLD_DIR ngsLD prune_graph THREADS_PER_CHROM FILTERED_PREFIX MAX_DIST MAX_KB_DIST R2

# Run the process in parallel for each chromosome, limiting to $MAX_JOBS
cat $CHROMOSOMES_FILE | parallel -j $MAX_JOBS process_chromosome {}

# Stitch all the unlinked positions into one file
cat $ngsLD_DIR/${FILTERED_PREFIX}_*_unlinked.pos > $ngsLD_DIR/${FILTERED_PREFIX}_unlinked.pos

echo "All chromosomes processed and unlinked SNPs combined."