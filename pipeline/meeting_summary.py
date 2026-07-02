import sys
import pickle
sys.path.append("/home/shatcher1/projects/PrimateAI-3D/pipeline")

from individual import Individual

with open("/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_individuals.pkl", "rb") as f:
    individuals = pickle.load(f)

# print summary for first 10 individuals
print(f"{'SAMPLE_ID':<12} {'PATHOGENIC':<12} {'BENIGN':<12} {'TOTAL VARIANTS'}")
print("-" * 50)
for ind in individuals[:10]:
    burden = ind.compute_burden()
    print(f"{ind.sample_id:<12} {burden['pathogenic']:<12} {burden['benign']:<12} {len(ind.variants)}")

# print dispersion across all individuals
print("\n--- Dispersion (chr22) ---")
dispersion = Individual.compute_dispersion(individuals)
print(f"Deleterious: {dispersion['pathogenic']:.4f}")
print(f"Benign:     {dispersion['benign']:.4f}")
print(f"\nTotal individuals: {len(individuals)}")