import os
import fastdfe as fd

species="atlantic_cod"

DATA_DIR = os.environ.get("DATA_DIR", "../../data")

fasta_path = f"{DATA_DIR}/{species}/genome/genome.fna"
gff_path = f"{DATA_DIR}/{species}/genome/annotation.gff"

n_ingroups=30
outgroups=["DRR622371", "SRR21531029"]

filtering_proportion=round(583684/872156,2)

vcf_dir = f"{DATA_DIR}/outgroups/7_polarize/{species}"
output_dir = f"{DATA_DIR}/outgroups/7_polarize/{species}"

# filepaths
input_vcf = os.path.join(vcf_dir, f"{species}_merge_outgroup_fastdfe.final.cds.vcf.gz")
output_vcf = os.path.join(output_dir, f"{species}_fastdfe.cds.aa.vcf.gz")

# chromosomes to process
chroms = [
    "NC_044048.1", "NC_044049.1", "NC_044050.1", "NC_044051.1",
    "NC_044052.1", "NC_044053.1", "NC_044054.1", "NC_044055.1",
    "NC_044056.1", "NC_044057.1", "NC_044058.1", "NC_044059.1",
    "NC_044060.1", "NC_044061.1", "NC_044062.1", "NC_044063.1",
    "NC_044064.1", "NC_044065.1", "NC_044066.1", "NC_044067.1",
    "NC_044068.1", "NC_044069.1", "NC_044070.1"
]

# Compute target sites for all chromosomes
target_sites = fd.Annotation.count_target_sites(gff_path, remove_overlaps=True, contigs=chroms)
n_target_sites = sum(target_sites.values()) * filtering_proportion
print(f"n_target_sites: {n_target_sites}")
# Create an adaptive prior that allows divergence (mono-allelic in the ingroup, divergence from outgroup)
prior = fd.KingmanPolarizationPrior(allow_divergence=True)
# annotator instance
ann = fd.Annotator(
    vcf=input_vcf,
    fasta=fasta_path,
    gff=gff_path,
    annotations=[fd.MaximumLikelihoodAncestralAnnotation(
        outgroups=outgroups,
        n_ingroups=n_ingroups,
        parallelize=True,
        n_target_sites=n_target_sites,
        max_sites=10000,
        prior=prior
    )],
output=output_vcf
)
# Run the annotation and plot polarization probs
ann.annotate()
ax = prior.plot(file=os.path.join(output_dir, f"polarization_plot.cds.png"), show=False)