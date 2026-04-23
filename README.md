# genomic-load-marine-fish

Code for the manuscript "Genomic load accumulation during parallel
colonization of the Baltic Sea", covering analyses of genetic load,
diversity and demographic history in five North Atlantic marine fish
species: Atlantic cod, Atlantic herring, European flounder, lumpfish
and turbot.

## Repository layout

```
scripts/
  Preprocessing/      Raw-read QC, trimming, mapping, deduplication, downsampling
  Variant_calling/    High-coverage (bcftools) and low+downsampled (ANGSD) calling,
                      depth distributions, ngsLD decay and pruning
  Outgroups/          Outgroup variant calling, ancestral-state probability
                      estimation (fastdfe), direction of selection (DoS),
                      R_xy with block jackknife, purifying-load estimation
  SLiM/               Forward simulations under the Kyriazis DFE framework
  psmc.sh	      Run PSMC to estimate demographic history trajectories for each species
data/                 Per-species inputs and intermediate files (not tracked)
figures/              Generated plots
```

## Data

Genomic data are not included in this repository. Raw sequencing reads will
be deposited under accessions listed in the manuscript.
