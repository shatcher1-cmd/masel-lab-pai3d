import sys
sys.path.append("/home/shatcher1/projects/PrimateAI-3D/pipeline")

import pickle
from individual import Individual

with open("/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_individuals.pkl", "rb") as f:
    individuals = pickle.load(f)

ind = individuals[0]
print(f"Sample ID: {ind.sample_id}")
print(f"Population: {ind.population}")

# also check how many individuals got population=None, which would mean a mismatch
none_count = sum(1 for i in individuals if i.population is None)
print(f"Individuals with no population assigned: {none_count}")

dispersion_gbr = Individual.compute_dispersion_by_population(individuals, "GBR")
print(dispersion_gbr)

# Memory Exploration

print("---------Memory Exploration---------")

total_variants = sum(len(ind.variants) for ind in individuals)
avg_variants = total_variants / len(individuals)

print(f"Total individuals: {len(individuals)}")
print(f"Total Variants objects across all individuals: {total_variants}")
print(f"Average variants per individual: {avg_variants:.1f}")

single_variant = individuals[0].variants[0]
print(f"Size of one Variants object: {sys.getsizeof(single_variant)} bytes")

# size of the list itself for one individual (just the pointers, not the objects)
single_individuals_list = individuals[0].variants
print(f"Size of one individual's variants list (pointers only): {sys.getsizeof(single_individuals_list)} bytes")

# size of the outer list of 3,202 Individual objects (pointers only)
print(f"Size of the individuals list itself (pointers only): {sys.getsizeof(individuals)} bytes")







