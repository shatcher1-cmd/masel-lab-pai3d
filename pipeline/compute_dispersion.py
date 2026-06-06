import csv
import numpy as np
import os

# paths
BURDEN_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_burden.tsv"
DISPERSION_OUTPUT_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_dispersion.tsv"

def compute_dispersion(burden_path, dispersion_output_path):
    pathogenic_counts = []
    benign_counts = []

    # read burden file and collect counts for all individuals
    with open(burden_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip header

        for row in reader:
            # column 1 is pathogenic burden, column 2 is benign burden
            pathogenic_counts.append(int(row[1]))
            benign_counts.append(int(row[2]))

    # convert to numpy arrays for vectorized math
    pathogenic_counts = np.array(pathogenic_counts)
    benign_counts = np.array(benign_counts)

    # compute mean and variance for pathogenic burden
    # ddof=1 means sample variance (divides by n-1)
    path_mean = np.mean(pathogenic_counts)
    path_variance = np.var(pathogenic_counts, ddof=1)
    path_dispersion = path_variance / path_mean

    # compute mean and variance for benign burden
    benign_mean = np.mean(benign_counts)
    benign_variance = np.var(benign_counts, ddof=1)
    benign_dispersion = benign_variance / benign_mean

    # print results to terminal
    print(f"--- Pathogenic ---")
    print(f"  Mean:       {path_mean:.4f}")
    print(f"  Variance:   {path_variance:.4f}")
    print(f"  Dispersion: {path_dispersion:.4f}")
    print(f"--- Benign ---")
    print(f"  Mean:       {benign_mean:.4f}")
    print(f"  Variance:   {benign_variance:.4f}")
    print(f"  Dispersion: {benign_dispersion:.4f}")

    # write results to TSV
    with open(dispersion_output_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["CATEGORY", "MEAN", "VARIANCE", "INDEX_OF_DISPERSION"])
        writer.writerow(["pathogenic", path_mean, path_variance, path_dispersion])
        writer.writerow(["benign", benign_mean, benign_variance, benign_dispersion])

    print(f"Dispersion results written to {dispersion_output_path}")

if __name__ == "__main__":
    compute_dispersion(BURDEN_PATH, DISPERSION_OUTPUT_PATH)