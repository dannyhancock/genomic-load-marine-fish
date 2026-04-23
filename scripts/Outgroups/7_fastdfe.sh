# This script is for merging ingroup and outgroup variants, using fastdfe to determine
# the ancestral state at each site, using snpEff to annotate the VCF and polarizing the VCF 
# to have the ancestral allele as reference

SPECIES="atlantic_cod"

# Focal VCF should already be hard filtered on QUAL and MQ and contain only biallelic SNPs on chroms
FOCAL_VCF_FULL=~/fish-genomics/data/${SPECIES}/High_coverage/${SPECIES}.QUAL30.MQ30.chroms.snps.bi.vcf.gz
PROTEIN_CODING_GENE_REGIONS=~/fish-genomics/data/${SPECIES}/genome/protein_coding_gene_regions.bed
CDS=~/fish-genomics/data/${SPECIES}/genome/cds_merged.bed
CHROMOSOMES=~/fish-genomics/data/${SPECIES}/High_coverage/chromosomes.txt
EXTRA_ARG="_fastdfe"

VARIANT_DIR=~/fish-genomics/data/outgroups/6_variant_calling/${SPECIES}
POLARIZE_DIR=~/fish-genomics/data/outgroups/7_polarize/${SPECIES}

# get mean depth of individuals for filtering
SPECIES_DEPTH=~/fish-genomics/data/${SPECIES}/High_coverage/$SPECIES.QUAL30.MQ30.chroms.snps.bi

vcftools --gzvcf $FOCAL_VCF_FULL \
--depth \
--out $SPECIES_DEPTH

#############
# CLOSEST OUTGROUP MUST BE OUTGROUP1, SECOND CLOSEST OUTGROUP2 AND SO ON!
#############
OUTGROUP1=DRR622371
OUTGROUP2=SRR21531029

# Extract Coding Sequences first
FOCAL_VCF=~/fish-genomics/data/${SPECIES}/High_coverage/${SPECIES}.QUAL30.MQ30.chroms.snps.bi.cds.vcf.gz
bcftools view $FOCAL_VCF_FULL -T $CDS -Oz -o $FOCAL_VCF --threads 4
bcftools index $FOCAL_VCF --threads 4

# concatenate ingroup variant VCF to ingroup invariant VCF and sort
### make sure INDELS are removed from both files first. ###
# We also only take coding sequences
INVARIANT_FULL=$VARIANT_DIR/focal_invar_with_${OUTGROUP1}_${OUTGROUP2}_var_intersect.vcf.gz
INVARIANT=$VARIANT_DIR/focal_invar_with_${OUTGROUP1}_${OUTGROUP2}_var_intersect.cds.vcf.gz

bcftools view $INVARIANT_FULL -T $CDS -Oz -o $INVARIANT --threads 4
bcftools index $INVARIANT --threads 4

# temp file for concatenation
tmpfocal=$(mktemp --suffix=.vcf.gz)

bcftools concat -a $FOCAL_VCF $INVARIANT -Ou \
	| bcftools sort -Oz -o $tmpfocal

# Count total sites (roughly) that are non-monomorphic across ingroup and outgroups for fastdfe
# This is used to determine the proportion of all bases in coding sequences to draw monomorphic sites from 
bcftools view $tmpfocal -H | wc -l > $POLARIZE_DIR/filtering/nbSNPs_merge_pre_filter.txt

###########################
## FILTERING and MERGING ##
###########################
mkdir -p $POLARIZE_DIR/filtering
# Filter ingroup VCF plus invariant sites for min and max depth, GQ and no missing data. Keeping only regions in protein coding genes.
maxDP=$(awk '{ sum += $3 } END { printf("%.0f", 2.5 * sum/(NR-1)) }' "$SPECIES_DEPTH.idepth")
minDP=5
minGQ=30
MISS_FILT=0.0

FILTERED_FOCAL_VCF="$POLARIZE_DIR/filtering/${SPECIES}.QUAL30.MQ30.chroms.snps.bi.minDP${minDP}.maxDP${maxDP}.minGQ${minGQ}.MISS${MISS_FILT}.cds.vcf.gz"

bcftools +setGT $tmpfocal \
 -- -t q -n . -i "FMT/DP < ${minDP} | FMT/GQ < ${minGQ}" \
 | bcftools +fill-tags \
 | bcftools view -T $CDS \
 -i "F_MISSING <= ${MISS_FILT} && MEAN(FMT/DP) <= ${maxDP}" \
 -m1 -M2 -Oz -o $FILTERED_FOCAL_VCF \
 --threads 4
 
bcftools index $FILTERED_FOCAL_VCF --threads 4
rm $tmpfocal

# Now filter outgroups
vcftools --gzvcf $VARIANT_DIR/${OUTGROUP1}.noindels.bi.vcf.gz \
--depth \
--out $VARIANT_DIR/${OUTGROUP1}.noindels.bi
vcftools --gzvcf $VARIANT_DIR/${OUTGROUP2}.noindels.bi.vcf.gz \
--depth \
--out $VARIANT_DIR/${OUTGROUP2}.noindels.bi

# Now filter outgroups
maxDP_OUTGROUP1=$(awk '{ sum += $3 } END { printf("%.0f", 2.5 * sum/(NR-1)) }' $VARIANT_DIR/${OUTGROUP1}.noindels.bi.idepth)
maxDP_OUTGROUP2=$(awk '{ sum += $3 } END { printf("%.0f", 2.5 * sum/(NR-1)) }' $VARIANT_DIR/${OUTGROUP2}.noindels.bi.idepth)

bcftools merge \
  $VARIANT_DIR/${OUTGROUP1}.noindels.bi.vcf.gz \
  $VARIANT_DIR/${OUTGROUP2}.noindels.bi.vcf.gz \
  -Oz -o $VARIANT_DIR/${OUTGROUP1}_${OUTGROUP2}_merged.noindels.bi.vcf.gz --threads 4

bcftools +setGT $VARIANT_DIR/${OUTGROUP1}_${OUTGROUP2}_merged.noindels.bi.vcf.gz \
 -- -t q -n . -i "FMT/GQ < ${minGQ} | FMT/DP[0] > ${maxDP_OUTGROUP1} | FMT/DP[1] > ${maxDP_OUTGROUP2}" \
 | bcftools +fill-tags \
 | bcftools view -T $CDS \
 -i "F_MISSING <= ${MISS_FILT}" \
 -m1 -M2 -Oz -o $POLARIZE_DIR/filtering/${OUTGROUP1}_${OUTGROUP2}_merged_filtered.cds.vcf.gz \
 --threads 4

bcftools index $POLARIZE_DIR/filtering/${OUTGROUP1}_${OUTGROUP2}_merged_filtered.cds.vcf.gz --threads 4

# then merge filtered outgroup VCF (including invariant sites in OUTGROUPS) to focal species
bcftools merge \
	$FILTERED_FOCAL_VCF \
	$POLARIZE_DIR/filtering/${OUTGROUP1}_${OUTGROUP2}_merged_filtered.cds.vcf.gz \
	-Oz -o $POLARIZE_DIR/filtering/${SPECIES}_merge_outgroup${EXTRA_ARG}.cds.vcf.gz \
	--threads 4

# Merging outgroups will result in lots of sites where ingroup is missing data. we 
# are not interested in these so we find sites with no missing data for ingroup
bcftools view -s "^${OUTGROUP1},${OUTGROUP2}" $POLARIZE_DIR/filtering/${SPECIES}_merge_outgroup${EXTRA_ARG}.cds.vcf.gz | \
  bcftools view -i "F_MISSING<=${MISS_FILT}" | \
  bcftools query -f '%CHROM\t%POS\n' > $POLARIZE_DIR/filtering/ingroup_nonmissing_sites${EXTRA_ARG}.cds.txt

# keep only sites with no missing data in the ingroup
bcftools view -T $POLARIZE_DIR/filtering/ingroup_nonmissing_sites${EXTRA_ARG}.cds.txt $POLARIZE_DIR/filtering/${SPECIES}_merge_outgroup${EXTRA_ARG}.cds.vcf.gz \
-m2 -M2 -Oz -o $POLARIZE_DIR/${SPECIES}_merge_outgroup${EXTRA_ARG}.final.cds.vcf.gz --threads 4

bcftools index $POLARIZE_DIR/${SPECIES}_merge_outgroup${EXTRA_ARG}.final.cds.vcf.gz --threads 4

# get missingness per indv to check completeness of outgroups
vcftools --gzvcf $POLARIZE_DIR/${SPECIES}_merge_outgroup${EXTRA_ARG}.final.cds.vcf.gz \
--missing-indv \
--out $POLARIZE_DIR/${SPECIES}_merge_outgroup${EXTRA_ARG}.cds

# output the number of snps after filtering
bcftools view $POLARIZE_DIR/${SPECIES}_merge_outgroup${EXTRA_ARG}.final.cds.vcf.gz -H | wc -l > $POLARIZE_DIR/filtering/nbSNPs_merge_post_filter${EXTRA_ARG}.txt

# use fastdfe to find and annotate ancestral states (AA and AA_prob).
python 7.2_annotate_ancestral.py

##########################
## ANNOTATE with snpEff ##
##########################
VCF_TO_ANNOTATE=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}.cds.aa.vcf.gz
ANNOTATED_VCF=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}.cds.ann.vcf.gz

# First, annotate the VCF with snpEff, so that the correct REF is used for annotation
SPECIES_DB="gadMor3.0"

java -jar ~/snpEff/snpEff.jar \
	-v $SPECIES_DB \
	$VCF_TO_ANNOTATE | bgzip > $ANNOTATED_VCF

bcftools index $ANNOTATED_VCF --threads 4

##################
## Polarization ##
##################

# keep sites where probability of the ancestral allele being called correctly is >= 60%
ANNOTATED_HIGH_PROB=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}.cds.aa0.6.ann.vcf.gz

bcftools view $ANNOTATED_VCF -i "AA_prob >= 0.6" \
-Oz -o $ANNOTATED_HIGH_PROB --threads 4

# Polarize the VCF to have the Ancestral Allele (AA) as the REF allele 
# and the Derived Allele as the Alternative Allele.
POLARIZED_VCF=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}_polarized.cds.ann.vcf.gz
bash ~/fish-genomics/scripts/outgroups/repolarize_vcf.sh $ANNOTATED_HIGH_PROB | bgzip -c > $POLARIZED_VCF

bcftools index $POLARIZED_VCF --threads 4
