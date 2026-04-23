#!/bin/bash
# repolarize_vcf.sh
# Author: Danny Hancock - daaha[at]dtu.dk
# Last updated: 15/04/2025
# This script repolarizes a VCF file so that the ancestral allele (from INFO/AA)
# becomes the reference. If the ancestral allele differs from REF, it swaps REF and ALT
# updates allele frequencies (AF) and updates genotype calls (swapping 0 and 1).
#
# Usage: bash repolarize_vcf.sh input.vcf.gz > output.vcf.gz

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 input.vcf" >&2
  exit 1
fi

input_vcf="$1"

zcat "$input_vcf" | awk 'BEGIN { FS="\t"; OFS="\t"; site_count=0; flip_count=0; missingAA=0 }
  # Print header lines unchanged.
  /^#/ { print; next }

  {
	site_count++;
	
    # Parse key fields:
    # $1: CHROM, $2: POS, $3: ID, $4: REF, $5: ALT, $8: INFO, $9+: sample fields.
    ref = $4;
    alt = $5;
    info = $8;
    ancestral = "";

    # Extract the ancestral allele from INFO field (assuming format AA=<allele>).
    if (match(info, /AA=([^;|]+)/, arr)) {
      ancestral = toupper(arr[1]);
    }
	
	# If the AA tag was not annotated “.”, skip this record entirely
    if (ancestral == ".") {
      missingAA++
	  next
    }

    # Otherwise if ancestral allele is missing or invalid, default to REF and warn.
    if (ancestral == "" || !(ancestral ~ /^[ACGT]$/)) {
      ancestral = ref;
      print "Warning: Invalid or missing ancestral allele at " $1 ":" $2 > "/dev/stderr";
    }

    # If the REF allele is already the ancestral allele, print the record unchanged.
    if (ancestral == ref) {
      print;
      next;
    }

    # Otherwise, swap REF and ALT.
    $4 = ancestral;
    $5 = ref;
	flip_count++;

    # update an AF field in INFO if present (swapping allele frequencies).
    if (match(info, /AF=([0-9\.]+)/, af_arr)) {
      old_af = af_arr[1] + 0;
      new_af = 1 - old_af;
      gsub("AF=" old_af, "AF=" new_af, info);
      $8 = info;
    }

    # For each sample genotype (columns 10 onward), flip the genotype calls.
    for (i = 10; i <= NF; i++) {
      split($i, fields, ":");
      gt = fields[1];
      sep = "";
      if (index(gt, "/") > 0)
        sep = "/";
      else if (index(gt, "|") > 0)
        sep = "|";
      
      if (sep != "") {
        split(gt, alleles, sep);
        for (j in alleles) {
          if (alleles[j] == "0")
            alleles[j] = "1";
          else if (alleles[j] == "1")
            alleles[j] = "0";
        }
        new_gt = alleles[1];
        for (k = 2; k in alleles; k++) {
          new_gt = new_gt sep alleles[k];
        }
        # Reattach the remaining genotype fields (if any).
        if (length(fields) > 1)
          new_gt = new_gt ":" substr($i, index($i, ":") + 1);
        $i = new_gt;
      }
    }
    print;
  }
  END { 
         if (site_count > 0) {
           printf("Number of flipped sites: %d out of %d (%.2f%%)\n", flip_count, site_count, (flip_count/site_count)*100) > "/dev/stderr"
		   printf("Number of records skipped due to AA==\".\": %d\n", missingAA) > "/dev/stderr"
         } else {
           print "No sites processed." > "/dev/stderr"
         }
  }'
  