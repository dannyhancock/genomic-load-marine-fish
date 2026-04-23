SPECIES="atlantic_cod"
FOLDER="fastdfe"
EXTRA_ARG="_fastdfe"

POLARIZE_DIR=~/fish-genomics/data/outgroups/7_polarize/${SPECIES}
DoS_DIR=~/fish-genomics/data/outgroups/8_DoS/${SPECIES}
LOAD_DIR=~/fish-genomics/data/outgroups/9_Load/${SPECIES}

INPUT_VCF=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}_polarized.cds.ann.vcf.gz

gene_regions=~/fish-genomics/data/${SPECIES}/genome/protein_coding_gene_regions.bed

IMPACT_LEVELS=("HIGH" "MODERATE" "LOW")
POPULATIONS=("baltic_sea" "kattegat" "north_sea")

NOFIXED=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}_polarized.nofixed.cds.ann.vcf.gz
OUTGROUP1=DRR622371
OUTGROUP2=SRR21531029
# Remove sites that are fixed derived across the entire species and therefore unlikely to be deleterious #
bcftools view $INPUT_VCF -s "^${OUTGROUP1},${OUTGROUP2}" | \
	bcftools +fill-tags -Ou -- -t AF | \
	bcftools view -i "AF != 1" -Oz -o $NOFIXED --threads 4
INPUT_VCF=$NOFIXED

# Select which genes we want to calculate load in 
direction="purifying"
gene_selection="all"
GENE_SELECT_DIR=$LOAD_DIR/${FOLDER}/${SPECIES}_${gene_selection}_${direction}
SAMPLE_SNPS_DIR=$GENE_SELECT_DIR/sample_snps
mkdir -p $SAMPLE_SNPS_DIR
OUTPUT_CSV=$GENE_SELECT_DIR/${SPECIES}_${gene_selection}_${direction}_nofixed_load.csv
echo "sample,pop,impact,regular,masked,realized,total" > "$OUTPUT_CSV"

for pop in "${POPULATIONS[@]}"; do
    echo "Processing population: $pop"
    
    # Population-specific files.
    pop_map="${LOAD_DIR}/inds_${pop}.txt"
    # purifying_genes="genes/${pop}${gene_selection}_pur.txt"
	
	# For the genes identified as purifying in the entire species
	purifying_genes="${LOAD_DIR}/${FOLDER}/genes/${SPECIES}_${gene_selection}_${direction}.txt"
    
    # Create a BED file containing only the purifying gene regions.
    awk 'NR==FNR { pur[$1] = 1; next }
    {
        split($4, fields, ";");
        gene = "";
        for(i in fields) {
            if(fields[i] ~ /^gene=/) {
                split(fields[i], kv, "=");
                gene = kv[2];
                break;
            }
        }
        if(gene in pur) print $0;
    }' "$purifying_genes" "$gene_regions" > "${GENE_SELECT_DIR}/${pop}_${gene_selection}_regions.bed"
    
    # Get list of samples for this population.
    SAMPLES=($(cat "$pop_map"))
    
    for sample in "${SAMPLES[@]}"; do
        echo "  Processing sample: $sample (Population: $pop)"
        sample_file="$SAMPLE_SNPS_DIR/${sample}_snps.txt.gz"
        
        # extract the SNP positions, GT and INFO/ANN field.
		bcftools view -T "${GENE_SELECT_DIR}/${pop}_${gene_selection}_regions.bed" -s "$sample" "$INPUT_VCF" | \
		bcftools query -f '%CHROM\t%POS\t[%GT]\t%INFO/ANN\n' | gzip > "$sample_file"
        
        # For each impact level, count SNPs based on the first annotation in the ANN field.
        for impact in "${IMPACT_LEVELS[@]}"; do
            echo "    Counting $impact impact variants..."
            counts=$(gunzip -c "$sample_file" | awk -v impact="$impact" '
            {
                # Fields: $1 = CHROM, $2 = POS, $3 = GT, $4 = INFO/ANN
                # Use only the first annotation (split on comma).
                split($4, annos, ",");
                split(annos[1], fields, "|");
                # Check if the third field (impact) exactly equals the target.
                if(fields[3] == impact) {
                    if($3 == "0/0") regular++;
                    else if($3 == "0/1" || $3 == "1/0") masked++;
                    else if($3 == "1/1") realized++;
                }
            }
            END {
                print (regular+0), (masked+0), (realized+0)
            }')
            
            # Read the counts into variables.
            read regular masked realized <<< "$counts"
            total=$(( regular + masked + realized ))
            
            # Write a CSV row for this sample, population, and impact level.
            echo "${sample},${pop},${impact},${regular},${masked},${realized},${total}" >> "$OUTPUT_CSV"
        done
    done
done

echo "High, moderate and low impact mutations in ${direction} genes counted for each population."

