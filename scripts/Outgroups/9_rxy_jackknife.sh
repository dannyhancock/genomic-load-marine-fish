jackknife_csv=~/fish-genomics/data/outgroups/9_Load/pairwise_rxy_jackknife.csv
echo "species,pop1,pop2,impact,Rxy,SE,Z,P" > "$jackknife_csv"

B=100   # number of blocks

for SPECIES in "atlantic_cod" "atlantic_herring" "european_flounder" "lumpfish" "turbot"; do
(
	echo $SPECIES
	if [ "$SPECIES" = "atlantic_cod" ]; then
	  populations=( "north_sea" "kattegat" "baltic_sea" )
	elif [ "$SPECIES" = "atlantic_herring" ]; then
	  populations=( "north_sea" "kattegat" "baltic_sea" )
	elif [ "$SPECIES" = "european_flounder" ]; then
	  populations=( "north_sea" "kattegat" "bornholm_basin" "baltic_sea" )
	elif [ "$SPECIES" = "lumpfish" ]; then
	  populations=( "skagerrak" "kattegat" "oresund" "baltic_sea" )
	elif [ "$SPECIES" = "turbot" ]; then
	  populations=( "north_sea" "kattegat" "baltic_sea" )
	else
	  echo "Unknown SPECIES: $SPECIES" >&2
	  exit 1
	fi
	
	echo "Processing $SPECIES"
	
	LOWCOV_DIR=~/fish-genomics/data/${SPECIES}/Low_coverage
	MIXED_DIR=~/fish-genomics/data/${SPECIES}/Mixed
	POLARIZED_VCF=~/fish-genomics/data/outgroups/7_polarize/${SPECIES}/${SPECIES}_fastdfe_polarized.nofixed.cds.ann.vcf.gz

	LOAD_DIR=~/fish-genomics/data/outgroups/9_Load/${SPECIES}
	rxy_DIR=${LOAD_DIR}/rxy

	mkdir -p ${rxy_DIR}

	MODEL="fastdfe"
	direction="purifying"
	gene_selection="all"
	GENE_SELECT_DIR=${LOAD_DIR}/${MODEL}/${SPECIES}_${gene_selection}_${direction} 
	

	# Extract CHR POS REF ALT and the SnpEff “IMPACT” field from the ANN tag:
	bcftools view $POLARIZED_VCF \
	-T ${GENE_SELECT_DIR}/${populations[0]}_${gene_selection}_regions.bed \
	| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' \
	| awk -F'\t' '{
		# $5 is the full ANN=... string; field[3] of that (pipe delimited) is IMPACT
		n = split($5,a,"|");
		imp = a[3];

		# print a “lookup key” of CHR:POS, then REF, ALT, IMPACT
		printf "%s:%s\t%s\t%s\t%s\n", $1, $2, $3, $4, imp
	  }' OFS='\t' \
	> ${rxy_DIR}/vcf_lookup.txt
	
	##############################################
    ### Per-population derived AF tables
    ##############################################
	
	for POP in "${populations[@]}"; do
	  echo "Running population: ${POP}"
	  MAFGZ="${MIXED_DIR}/${POP}_lc_and_downsampled.cds.masked.mafs.gz"
	  OUTTAB="${rxy_DIR}/${POP}_lc.tab"

	  max_nInd=$(zcat "${MAFGZ}" | awk -F'\t' 'NR>1 && $7~/^[0-9]/ && $7>m{m=$7} END{print m}')
	  echo "max nInd=${max_nInd}"

	  zcat "${MAFGZ}" \
		| awk -v -F'\t' '
		  NR>1 {
			key = $1 ":" $2
			print key, $3, $4, $5, $6
		  }
		' OFS='\t' \
	  > "${rxy_DIR}/${POP}_lc.freqs.txt"

	  # Join to vcf_lookup.txt on “key” to grab (anc,der,impact) per site
	  join -t $'\t' \
		   <( sort ${rxy_DIR}/${POP}_lc.freqs.txt ) \
		   <( sort ${rxy_DIR}/vcf_lookup.txt ) \
	  | awk -F'\t' '{
		  # fields after join:
		  #   $1 = key (CHR:POS)
		  #   $2 = major
		  #   $3 = minor
		  #   $4 = ref
		  #   $5 = knownEM    (this is P(minor allele))
		  #	  $6 = anc
		  #   $7 = der
		  #   $8 = impact

		  key    = $1
		  major  = $2
		  minor  = $3
		  ref    = $4
		  em     = $5
		  anc    = $6
		  der    = $7
		  imp    = $8

		  # Compute derived allele freq (ALT=derived in the polarized VCF)
		  if (minor == der) {
			daf = em
		  }
		  else if (minor == anc) {
			daf = 1 - em
		  }
		  else {
			# The “minor” allele in MAF didn’t match either ANC or DER → skip it
			next
		  }

		  # split key into chr, pos
		  split(key, K, ":")
		  chr = K[1]
		  pos = K[2]

		  # Print exactly the same 5 column “.tab” as before:
		  #   key, CHR, POS, IMPACT, derivedAF
		  print key, chr, pos, imp, daf
		}' OFS='\t' \
	  > ${OUTTAB}

	  # Cleanup intermediate:
	  rm ${rxy_DIR}/${POP}_lc.freqs.txt
	done

	##############################################
	### Pairwise Rxy + JACKKNIFE (single-output)
	##############################################

	pair_dir="${rxy_DIR}/pairwise"
	mkdir -p "${pair_dir}"

	N=${#populations[@]}
	for (( i=0; i<N; i++ )); do
		popA=${populations[i]}
		for (( j=i+1; j<N; j++ )); do
			popB=${populations[j]}

			echo "=== Rxy for ${popA} vs ${popB} ==="

			merged="${pair_dir}/${popA}_${popB}_merged.tab"

			# SNP-level merge
			join -t $'\t' -1 1 -2 1 \
				 -o '1.2,1.3,1.4,1.5,2.5' \
				 <(sort "${rxy_DIR}/${popA}_lc.tab") \
				 <(sort "${rxy_DIR}/${popB}_lc.tab") \
			  > "${merged}"

			# Re-sort numerically by CHR then POS for contiguous block assignment
			sort -t$'\t' -k1,1V -k2,2n "${merged}" -o "${merged}"

			##################################################
			## FULL-DATA Rxy (needed for Z and P)
			##################################################
			echo "Calculating full Rxy statistic on entire dataset"
			rxy_out="${pair_dir}/${popA}_${popB}_full.rxy"

			awk -F'\t' '
			  {
				imp=$3; f1=$4; f2=$5;
				c1[imp]+=f1*(1-f2);
				c2[imp]+=f2*(1-f1);
			  }
			  END {
				for (i in c1)
				  printf "%s\t%.10f\t%.10f\t%.10f\n", i, c1[i], c2[i], c1[i]/c2[i];
			  }
			' "${merged}" > "${rxy_out}"

			# load Rfull into shell variables
			declare -A Rfull

			while IFS=$'\t' read -r IMP c1 c2 R; do
				Rfull[$IMP]=$R
			done < "${rxy_out}"

			##################################################
			### Assign SNPs to jackknife blocks
			##################################################
			echo "Assigning SNPs to blocks"
			merged_blocks="${merged}.blocks"
			# First get number of lines
			total=$(wc -l < "${merged}")

			# Now assign blocks
			awk -v B=$B -v total="$total" -F'\t' '
			{
				blk = int( (NR-1) * B / total ) + 1;
				print $0, blk;
			}' OFS='\t' "${merged}" > "${merged_blocks}"

			# Count SNPs per impact class per block (for weighted jackknife)
			block_counts="${pair_dir}/${popA}_${popB}_block_counts.txt"
			awk -F'\t' '{
				imp = $3; blk = $6
				count[blk, imp]++
				total[imp]++
			} END {
				for (key in count) {
					split(key, a, SUBSEP)
					printf "%s\t%s\t%d\t%d\n", a[1], a[2], count[key], total[a[2]]
				}
			}' "${merged_blocks}" > "${block_counts}"

			##################################################
			### JACKKNIFE LOOP
			##################################################
			echo "Starting Jacknife"
			
			jk_file="${pair_dir}/${popA}_${popB}_jackknife_reps_10.txt"
			> "$jk_file"

			for b in $(seq 1 $B); do
				echo "Leaving out block $b"
				awk -v b=$b -F'\t' '$6 != b' "${merged_blocks}" > "${pair_dir}/${popA}_${popB}_tmp.no_block"

				awk -F'\t' -v blk=$b '
					{
						imp=$3; f1=$4; f2=$5;
						c1[imp]+=f1*(1-f2);
						c2[imp]+=f2*(1-f1);
					}
					END {
						for (i in c1)
							printf "%d\t%s\t%.10f\t%.10f\t%.10f\n",
								   blk, i, c1[i], c2[i], c1[i]/c2[i];
					}
				' "${pair_dir}/${popA}_${popB}_tmp.no_block" >> "$jk_file"
			done

			rm "${pair_dir}/${popA}_${popB}_tmp.no_block"

			##################################################
			### COMPUTE SE, Z, P FOR EACH IMPACT CLASS
			##################################################

			for IMP in HIGH MODERATE LOW; do

				R=${Rfull[$IMP]}
				if [[ "$R" == "NA" ]]; then
					echo "$SPECIES,$popA,$popB,$IMP,NA,NA,NA,NA" >> "${jackknife_csv}.${SPECIES}"
					continue
				fi

				SE=$(awk -v IMP="$IMP" -v B=$B -v Rfull="$R" -F'\t' '
					# File 1: block_counts (block, impact, n_b, n_total)
					FILENAME == ARGV[1] && $2 == IMP {
						n_b[$1] = $3
						n_total = $4
					}
					# File 2: jackknife_reps (block, impact, c1, c2, Rxy)
					FILENAME == ARGV[2] && $2 == IMP {
						r_loo[$1] = $5
					}
					END {
						if (n_total == 0) { print "NA"; exit }

						# Compute pseudovalues
						nb = 0
						for (b = 1; b <= B; b++) {
							if (!(b in n_b) || n_b[b] == 0) continue
							h = n_total / n_b[b]
							pv[b] = h * Rfull - (h - 1) * r_loo[b]
							h_arr[b] = h
							nb++
						}

						if (nb < 2) { print "NA"; exit }

						# Mean pseudovalue
						sum = 0
						for (b in pv) sum += pv[b]
						pv_mean = sum / nb

						# Weighted variance
						var = 0
						for (b in pv) {
							var += ((pv[b] - pv_mean)^2) / (h_arr[b] - 1)
						}
						var /= nb

						print sqrt(var)
					}
				' "${block_counts}" "$jk_file")

				if [ "$SE" = "NA" ]; then
					echo "$SPECIES,$popA,$popB,$IMP,$R,NA,NA,NA" >> "${jackknife_csv}.${SPECIES}"
					continue
				fi

				Z=$(awk -v R=$R -v SE=$SE 'BEGIN{print (R-1)/SE}')
				# P-value left blank for calculation in R
				P=""

				echo "$SPECIES,$popA,$popB,$IMP,$R,$SE,$Z,$P" >> "${jackknife_csv}.${SPECIES}"

			done

			echo "--> Finished ${popA} vs ${popB}"

		done
	done

) &
done

# Wait for all species to finish
wait

# Concatenate per-species results into final CSV
for SPECIES in "atlantic_cod" "atlantic_herring" "european_flounder" "lumpfish" "turbot"; do
	if [ -f "${jackknife_csv}.${SPECIES}" ]; then
		cat "${jackknife_csv}.${SPECIES}" >> "$jackknife_csv"
		rm "${jackknife_csv}.${SPECIES}"
	fi
done

echo "All species complete."
