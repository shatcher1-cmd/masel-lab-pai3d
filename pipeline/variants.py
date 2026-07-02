class Variants:
    def __init__(self, chromosome, gene, position, reference_allele, alternate_allele, genotype_copies, pai3d_score):
        self.chromosome = chromosome
        self.gene = gene
        self.position = position
        self.reference_allele = reference_allele
        self.alternate_allele = alternate_allele
        self.genotype_copies = genotype_copies
        self.pai3d_score = pai3d_score
        self.pai3d_prediction = "pathogenic" if pai3d_score > 0.821 else "benign"