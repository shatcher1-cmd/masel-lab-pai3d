import csv
import os

# paths
ANNOTATED_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_annotated.tsv"
BURDEN_OUTPUT_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_burden.tsv"

def parse_genotype(genotype):
    # replace | with / so we handle both phased (0|1) and unphased (0/1) formats
    # then split on / to get individual alleles -> ["0", "1"]
    alleles = genotype.replace("|", "/").split("/")
    # count how many alleles are "1" (the alternate/mutant allele)
    # 0|0 -> 0, 0|1 or 1|0 -> 1, 1|1 -> 2
    return sum(1 for a in alleles if a == "1")

def compute_burden(annotated_path, burden_output_path):
    with open(annotated_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        
        # read header row to extract sample IDs
        # header = ["CHROM", "POS", "REF", "ALT", "PREDICTION", "HG00096", "HG00097", ...]
        # sample IDs start at index 5
        header = next(reader)
        sample_ids = header[5:]

        # initialize burden counters for every individual at 0
        # these are dictionaries: {sample_id: count}
        pathogenic_burden = {sample: 0 for sample in sample_ids}
        benign_burden = {sample: 0 for sample in sample_ids}

        for row in reader:
            # column 4 is the PAI3D prediction: "benign" or "pathogenic"
            prediction = row[4]
            # columns 5 onwards are the genotype for each individual
            genotypes = row[5:]

            # zip pairs each sample ID with its genotype for this variant
            # e.g. ("HG00096", "0|1"), ("HG00097", "0|0"), ...
            for sample, genotype in zip(sample_ids, genotypes):
                # get copy count (0, 1, or 2) for this individual at this variant
                count = parse_genotype(genotype)
                # add to the correct burden counter based on prediction
                if prediction == "pathogenic":
                    pathogenic_burden[sample] += count
                elif prediction == "benign":
                    benign_burden[sample] += count

    # write output TSV with one row per individual
    with open(burden_output_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["SAMPLE_ID", "PATHOGENIC_BURDEN", "BENIGN_BURDEN"])
        for sample in sample_ids:
            writer.writerow([sample, pathogenic_burden[sample], benign_burden[sample]])

    print(f"Burden counts written for {len(sample_ids)} individuals")

if __name__ == "__main__":
    compute_burden(ANNOTATED_PATH, BURDEN_OUTPUT_PATH)
    print(f"Output written to {BURDEN_OUTPUT_PATH}")