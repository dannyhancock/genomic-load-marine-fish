############################
## Direction of Selection ##
############################
SPECIES="atlantic_cod"

POLARIZE_DIR=~/fish-genomics/data/outgroups/7_polarize/${SPECIES}
DoS_DIR=~/fish-genomics/data/outgroups/8_DoS/${SPECIES}

CDS=~/fish-genomics/data/${SPECIES}/genome/cds_merged.bed
CHROMOSOMES=~/fish-genomics/data/${SPECIES}/High_coverage/chromosomes.txt
EXTRA_ARG="_fastdfe"

FOLDER="fastdfe"
mkdir -p $DoS_DIR/$FOLDER
ANN_VCF=$POLARIZE_DIR/${SPECIES}${EXTRA_ARG}_polarized.cds.ann.vcf.gz

# Now we need to calculate the Direction of Selection for each gene, separately for each population
populations=(baltic_sea kattegat north_sea)

for POP in "${populations[@]}"; do
    echo "Processing population: $POP"
    # First, output the gene info with newly calculated allele frequencies.
    bcftools view -T "$CDS" -S $DoS_DIR/inds_"$POP".txt "$ANN_VCF" | \
        bcftools +fill-tags -Ou -- -t AF | \
        bcftools query -f '%CHROM\t%POS\t%INFO/ANN\t%INFO/AF\n' | \
     # Parse the gene table to extract the ANN field.
    awk -F "\t" '{
        # Split the ANN field on commas; take the first annotation.
        split($3, ann_arr, ",");
        split(ann_arr[1], fields, "|");
        gene = fields[5];        # Gene_ID
        var_type = fields[2];      # Variant type (e.g., may contain "synonymous_variant" or "missense_variant")
        impact = fields[3];        # Impact (e.g., HIGH, MODERATE, LOW, MODIFIER)
        print $1, $2, gene, var_type, impact, $4;
    }' OFS="\t" > "$DoS_DIR/$FOLDER/${POP}_genes_parsed.txt"

    # Calculate DoS per gene.
    awk -F "\t" '
    {
        gene = $3;
        all_genes[gene] = 1;
        af = $6 + 0;  # force numeric conversion of AF
        # Fixed differences (AF == 1)
        if (af == 1) {
			if ($4 ~ /missense_variant/)
				dn[gene]++;
			else if ($4 ~ /synonymous_variant/)
				ds[gene]++;
		}
        # Polymorphisms (0 < AF < 1)
        else if (af > 0 && af < 1) {
            if ($4 ~ /missense_variant/)
                pn[gene]++;
            else if ($4 ~ /synonymous_variant/)
                ps[gene]++;
        }
    }
    END {
        print "gene", "Dn", "Ds", "Pn", "Ps", "DoS";
        for (g in all_genes) {
            d_n = (dn[g] ? dn[g] : 0);
            d_s = (ds[g] ? ds[g] : 0);
            p_n = (pn[g] ? pn[g] : 0);
            p_s = (ps[g] ? ps[g] : 0);
            # if ((p_n + d_n) >= 4) {
                denom1 = d_n + d_s;
                denom2 = p_n + p_s;
                sub_rate = (denom1 > 0 ? d_n/denom1 : 0);
                poly_rate = (denom2 > 0 ? p_n/denom2 : 0);
                dos = sub_rate - poly_rate;
                print g, d_n, d_s, p_n, p_s, dos;
            # }
        }
    }' OFS="\t" "$DoS_DIR/$FOLDER/${POP}_genes_parsed.txt" > "$DoS_DIR/$FOLDER/${POP}_all_dos.txt"
done

## FULL SPECIES - just remove outgroups ##
POP="atlantic_cod"
echo "Processing population: $POP"
bcftools view -T $CDS -s "^${OUTGROUP1},${OUTGROUP2}" "$ANN_VCF" | \
        bcftools +fill-tags -Ou -- -t AF | \
        bcftools query -f '%CHROM\t%POS\t%INFO/ANN\t%INFO/AF\n' | \
# Parse the gene table to extract the ANN field.
awk -F "\t" '{
	# Split the ANN field on commas; take the first annotation.
	split($3, ann_arr, ",");
	split(ann_arr[1], fields, "|");
	gene = fields[5];        # Gene_ID
	var_type = fields[2];      # Variant type (e.g., may contain "synonymous_variant" or "missense_variant")
	impact = fields[3];        # Impact (e.g., HIGH, MODERATE, LOW, MODIFIER)
	print $1, $2, gene, var_type, impact, $4;
}' OFS="\t" > "$DoS_DIR/$FOLDER/${POP}_genes_parsed.txt"

# Calculate DoS per gene.
# for the entire species, we can look for sites that are all 0/0 or all 1/1 as these are 
# fixed differences (the VCF contains variant sites, so if the ingroup are fixed 0/0 the 
# outgroup must be variant
awk -F "\t" '
{
	gene = $3;
	all_genes[gene] = 1;
	af = $6 + 0;  # force numeric conversion of AF
	# Fixed differences (AF == 1)
	if (af == 1 || af == 0) {
		if ($4 ~ /missense_variant/)
			dn[gene]++;
		else if ($4 ~ /synonymous_variant/)
			ds[gene]++;
	}
	# Polymorphisms (0 < AF < 1)
	else if (af > 0 && af < 1) {
		if ($4 ~ /missense_variant/)
			pn[gene]++;
		else if ($4 ~ /synonymous_variant/)
			ps[gene]++;
	}
}
END {
	print "gene", "Dn", "Ds", "Pn", "Ps", "DoS";
	for (g in all_genes) {
		d_n = (dn[g] ? dn[g] : 0);
		d_s = (ds[g] ? ds[g] : 0);
		p_n = (pn[g] ? pn[g] : 0);
		p_s = (ps[g] ? ps[g] : 0);
		# if ((p_n + d_n) >= 4) {
			denom1 = d_n + d_s;
			denom2 = p_n + p_s;
			sub_rate = (denom1 > 0 ? d_n/denom1 : 0);
			poly_rate = (denom2 > 0 ? p_n/denom2 : 0);
			dos = sub_rate - poly_rate;
			print g, d_n, d_s, p_n, p_s, dos;
		# }
	}
}' OFS="\t" "$DoS_DIR/$FOLDER/${POP}_genes_parsed.txt" > "$DoS_DIR/$FOLDER/${POP}_all_dos_0s.txt"

# Standard MKtest
for POP in "${populations[@]}"; do
    echo "Processing population: $POP"
    # # First output the gene info with allele frequencies
	
	# STEP PROBABLY ALREADY DONE ABOVE
    # bcftools view -R "$GENE_REGIONS" -S inds_"$POP".txt "$ANN_VCF" | \
        # bcftools +fill-tags -Ou -- -t AF | \
        # bcftools query -f '%CHROM\t%POS\t%INFO/ANN\t%INFO/AF\n' -o "${FOLDER}/${MODEL}/${POP}_gene_table.txt"

    awk -F "\t" '
    {
        gene = $3;
        all_genes[gene] = 1;
        af = $6 + 0;
        # For fixed differences (AF == 1):
        if (af == 1) {
            if ($4 ~ /missense_variant/)
                dn[gene]++;
            else if ($4 ~ /synonymous_variant/)
                ds[gene]++;
        }
        # For polymorphisms (0 < AF < 1):
        else if (af > 0 && af < 1) {
            if ($4 ~ /missense_variant/)
                pn[gene]++;
            else if ($4 ~ /synonymous_variant/)
                ps[gene]++;
        }
    }
    END {
        print "gene", "Dn", "Ds", "Pn", "Ps", "NI";
        for (g in all_genes) {
            d_n = (dn[g] ? dn[g] : 0);
            d_s = (ds[g] ? ds[g] : 0);
            p_n = (pn[g] ? pn[g] : 0);
            p_s = (ps[g] ? ps[g] : 0);
            if ((p_n + d_n) >= 4) {
                if (p_s > 0 && d_n > 0)
                    NI = (p_n * d_s) / (p_s * d_n);
                else
                    NI = "NA";
                print g, d_n, d_s, p_n, p_s, NI;
            }
        }
    }' OFS="\t" "${FOLDER}/${MODEL}/${POP}_genes_parsed.txt" > "${FOLDER}/${MODEL}/${POP}_NI.txt"
done
