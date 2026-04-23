#!/bin/sh
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N ngsLD_decay_atlantic_cod
#PBS -e /home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/2_ngsLD/logs/ngsLD_decay_atlantic_cod.err
#PBS -o /home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/2_ngsLD/logs/ngsLD_decay_atlantic_cod.out
#PBS -m n
#PBS -l nodes=1:ppn=20
#PBS -l mem=75gb
#PBS -l walltime=02:00:00

set -e

module purge
module load tools
module load bcftools/1.20
module load gcc/12.2.0
module load intel/basekit/INITIALIZE/2023.0.0
module load intel/basekit/mkl/2023.0.0
module load R/4.4.0
module load gsl/2.7

##############
## Directories
##############
species="atlantic_cod"
ANGSD_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/1_ANGSD/atlantic_cod"
ngsLD_DIR="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/2_ngsLD/atlantic_cod"
mkdir --parents $ngsLD_DIR

CHROMS=$ANGSD_DIR/chromosomes.txt

# # Get min and max DP from distributions
depth_limits="/home/projects/dp_00007/data/BlueBioClimate/Data/Mixed/0_depth/depth_limits.txt"
minDepth=$(awk -v species="$species" '$1 == species {print $2}' "$depth_limits")
maxDepth=$(awk -v species="$species" '$1 == species {print $3}' "$depth_limits")

###########
## Software
###########
ngsLD="/home/projects/dp_00007/data/BlueBioClimate/Data/Software/ngsLD"
BCF_pref="$ANGSD_DIR/atlantic_cod_lc_and_downsampled_realigned.minDP${minDepth}.maxDP${maxDepth}.minInd${MIN_IND}"

# Check if the VCF is indexed with a CSI file
if [ ! -f "${BCF_pref}.bcf.csi" ]; then
    echo "CSI index not found. Indexing the VCF file..."
    bcftools index -c ${BCF_pref}.bcf --threads 20
else
    echo "CSI index already exists."
fi

######################
# Extract snps only on linkage groups, removing unmapped scaffolds/mtdna
######################
# Keep only chroms/linkage groups
bcftools view ${BCF_pref}.bcf -r $(paste -sd, $CHROMS) --threads 20 -Oz -o ${BCF_pref}.chroms.vcf.gz

# Repeat for Beagle file
chrom_list=$(paste -sd '|' $CHROMS)
zcat ${BCF_pref}.beagle.gz | awk -v chrom_list="$chrom_list" 'NR==1 || $1 ~ "^("chrom_list")_"' | gzip > ${BCF_pref}.chroms.beagle.gz

# Repeat for MAFS file
zcat ${BCF_pref}.mafs.gz | awk -v chrom_list="$chrom_list" 'NR==1 || $1 ~ "^("chrom_list")$"' | gzip > ${BCF_pref}.chroms.mafs.gz

#####################
## Prepare stripped BEAGLE and POS files for ngsLD
#####################
# Prepare a genotype likelihood file by removing the header row and the first three columns (i.e. positions, major allele, minor allele).
zcat ${BCF_pref}.chroms.beagle.gz | tail -n +2 | cut -f 4- | gzip > ${BCF_pref}.chroms.stripped.beagle.gz

# Prepare a pos file by removing the header and taking just the first two columns.
zcat ${BCF_pref}.chroms.mafs.gz | tail -n +2 | cut -f 1,2 | sed 's/:/_/g'| gzip > ${BCF_pref}.chroms.stripped.pos.gz

# extract no. individuals and sites for ngsLD
zcat ${BCF_pref}.chroms.stripped.beagle.gz | head -n 1 | awk '{print NF/3}' > $ANGSD_DIR/no_individuals.txt
zcat ${BCF_pref}.chroms.stripped.beagle.gz | wc -l > $ANGSD_DIR/no_sites.txt

n_inds=$(cat "$ANGSD_DIR/no_individuals.txt")
n_sites=$(cat "$ANGSD_DIR/no_sites.txt")

# #######################
# ## run sampled ngsLD ##
# #######################
$ngsLD/ngsLD \
--geno ${BCF_pref}.chroms.stripped.beagle.gz \
--pos ${BCF_pref}.chroms.stripped.pos.gz \
--probs \
--n_ind $n_inds \
--n_sites $n_sites \
--max_kb_dist 20 \
--min_maf 0.05 \
--n_threads 20 \
--rnd_sample 0.05  \
--out $ngsLD_DIR/atlantic_cod_samp_20k.ld

echo $ngsLD_DIR/atlantic_cod_samp_20k.ld > $ngsLD_DIR/ld_file_list_20k.txt

# run LD_decay script
Rscript --vanilla --slave $ngsLD/scripts/fit_LDdecay.R --ld_files $ngsLD_DIR/ld_file_list_20k.txt --out $ngsLD_DIR/atlantic_cod_samp_ld_plot_20k.pdf --ld r2 --fit_level 3

# delete temp stripped files
rm ${BCF_pref}.chroms.stripped.beagle.gz ${BCF_pref}.chroms.stripped.pos.gz

echo "LD Decay estimated"
