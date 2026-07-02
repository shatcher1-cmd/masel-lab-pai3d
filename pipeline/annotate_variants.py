import gzip
import pickle
import sys
import os

sys.path.append("/home/shatcher1/projects/PrimateAI-3D/pipeline")

from build_dict import build_dict
from variants import Variants
from individual import Individual

# paths
PAI3D_CHR22_PATH = "/home/shatcher1/projects/PrimateAI-3D/data/PrimateAI-3D.hg38.chr22.txt.gz"
PICKLE_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_variant_dict.pkl"
VCF_PATH = "/home/shatcher1/projects/PrimateAI-3D/data/1000g/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
OUTPUT_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_individuals.pkl"
POPULATION_PATH = "/home/shatcher1/projects/PrimateAI-3D/data/1000g/20130606_g1k_3202_samples_ped_population.txt"

def load_population_map(ped_path):
    # reads the 1000G ped/population file and builds a dict mapping sample_id -> population
    population_map = {}  # Initializes dictionary
    with open(ped_path, "r") as f:
        next(f)  # skip header row
        for line in f:
            fields = line.strip().split() # space delimited (not tab delimited)
            sample_id = fields[1]     # SampleID column
            population = fields[5]    # Population column
            population_map[sample_id] = population
    return population_map

# ==============================================================================
# STRAND ALIGNMENT / ALLELE MATCHING LOGIC
# PAI3D predictions use genomic forward-strand alleles (non_flipped_ref/_alt).
# However, the 1000 Genomes VCF may record variants on the reverse strand.
#
# To avoid falsely skipping valid matches due to string mismatching, the
# COMPLEMENT mapping allows the pipeline to check for reverse-complement matches
# on the fly.
# ==============================================================================

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def parse_genotype(genotype):
    # replace | with / so we handle both phased (0|1) and unphased (0/1) formats
    # then split on / to get individual alleles -> ["0", "1"]
    alleles = genotype.replace("|", "/").split("/")
    # count how many alleles are "1" (the alternate/mutant allele)
    # 0|0 -> 0, 0|1 or 1|0 -> 1, 1|1 -> 2
    return sum(1 for a in alleles if a == "1")

def annotate_variants(vcf_path, output_path, variant_dict, population_map):

    individuals = {}  # maps sample_id -> Individual object

    with gzip.open(vcf_path, "rt") as vcf:

        found = 0
        skipped = 0

        for line in vcf:

            # skip metadata lines
            if line.startswith("##"):
                continue

            # header line -> extract sample IDs and build one Individual per sample
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                sample_ids = fields[9:]

                # create one Individual object per sample, stored by sample_id
                for sample_id in sample_ids:
                    population = population_map.get(sample_id) # None if sample not found in file
                    individuals[sample_id] = Individual(sample_id, population = population)
                continue

            # data lines -> look up variant in PAI3D dict
            fields = line.strip().split("\t")
            key = (fields[0], fields[1], fields[3], fields[4])

            lookup = variant_dict.get(key)

            if lookup is None:
                # try complement strand
                ref_comp = COMPLEMENT.get(fields[3])
                alt_comp = COMPLEMENT.get(fields[4])

                if ref_comp and alt_comp:
                    flipped_key = (fields[0], fields[1], ref_comp, alt_comp)
                    lookup = variant_dict.get(flipped_key)

            if lookup is None:
                skipped += 1
                continue

            # unpack the tuple stored in the dict by build_dict
            gene_name, pai3d_score, pai3d_prediction = lookup

            # loop over each individual's genotype for this variant
            genotypes = fields[9:]
            for sample_id, genotype in zip(sample_ids, genotypes):
                copies = parse_genotype(genotype)

                # only create a VariantRecord if the individual carries at least one copy
                if copies > 0:
                    record = Variants(
                        chromosome=fields[0],
                        gene=gene_name,
                        position=int(fields[1]),
                        reference_allele=fields[3],
                        alternate_allele=fields[4],
                        genotype_copies=copies,
                        pai3d_score=pai3d_score
                    )
                    individuals[sample_id].add_variant(record)

            found += 1

    print(f"Variants matched: {found}")
    print(f"Variants skipped: {skipped}")

    # convert dict to list of Individual objects and pickle
    individual_list = list(individuals.values())

    with open(output_path, "wb") as f:
        pickle.dump(individual_list, f)

    print(f"Pickled {len(individual_list)} individuals to {output_path}")
    return individual_list


if __name__ == "__main__":
    variant_dict = build_dict(PAI3D_CHR22_PATH, PICKLE_PATH)
    population_map = load_population_map(POPULATION_PATH)
    annotate_variants(VCF_PATH, OUTPUT_PATH, variant_dict, population_map)