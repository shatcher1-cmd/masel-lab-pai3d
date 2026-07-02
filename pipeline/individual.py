import numpy as np
from variants import Variants

class Individual:
    def __init__(self, sample_id, population):
        self.sample_id = sample_id
        self.population = population
        self.variants = []          # List of Variants objects, populated by add_variant()

    def add_variant(self, variant):
        # appends a single Variants object to this individual's variant list
        self.variants.append(variant)

    def compute_burden(self, threshold=0.821):
        # counts total alternate allele copies across all variants scored as pathogenic
        # genotype_copies is 0, 1, or 2 - so this sums actual dosage, not just presence
        pathogenic_copies = 0
        benign_copies = 0
        for variant in self.variants:
            if variant.pai3d_score > threshold:
                pathogenic_copies += variant.genotype_copies
            else:
                benign_copies += variant.genotype_copies
        return {
            "pathogenic": pathogenic_copies,
            "benign": benign_copies
        }

    def compute_dispersion(individuals, threshold=0.821):
        # computes variance/mean ratio across a list of Individual objects
        # called on a list, not a single individual
        pathogenic_burdens = []
        benign_burdens = []
        # collect each individuals burden counts into lists
        for individual in individuals:
            burden = individual.compute_burden(threshold)
            pathogenic_burdens.append(burden["pathogenic"])
            benign_burdens.append(burden["benign"])
        # convert to numpy arrays for variance and mean calculations
        pathogenic_burdens = np.array(pathogenic_burdens)
        benign_burdens = np.array(benign_burdens)
        # variance / mean = index of dispersion
        # ddof=1 means sample variance (divides by n-1)
        pathogenic_dispersion = np.var(pathogenic_burdens, ddof=1) / np.mean(pathogenic_burdens)
        benign_dispersion = np.var(benign_burdens, ddof=1) / np.mean(benign_burdens)
        return {
            "pathogenic": pathogenic_dispersion,
            "benign": benign_dispersion
        }

    def compute_dispersion_by_population(individuals, population, threshold=0.821):
        # filters to only individuals matching the given population label,
        # then runs the same dispersion calculation as compute_dispersion
        filtered = []
        for ind in individuals:
            if ind.population == population:
                filtered.append(ind)
        if len(filtered) == 0:
            return None
        return Individual.compute_dispersion(filtered, threshold)