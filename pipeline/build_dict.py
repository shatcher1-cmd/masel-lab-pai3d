import gzip
import pickle
import os
import sys

# base directories for HPC
PAI3D_DIR = "/groups/masel/spencer_pai3d/data/pai3d"
PICKLE_DIR = "/groups/masel/spencer_pai3d/output"

def build_dict(chr_num):

    pai3d_path = os.path.join(PAI3D_DIR, f"PrimateAI-3D.hg38.chr{chr_num}.txt.gz")
    pickle_path = os.path.join(PICKLE_DIR, f"chr{chr_num}_variant_dict.pkl")

    # load and return the cached dictionary if the pickle file already exists
    if os.path.exists(pickle_path):
        print(f"Loading existing dictionary from pickle: {pickle_path}")
        with open(pickle_path, "rb") as f:
            return pickle.load(f)

    print(f"Building dictionary from PAI3D file: {pai3d_path}")

    variant_dict = {}

    with gzip.open(pai3d_path, "rt") as f:
        next(f)  # skip header
        for line in f:
            fields = line.strip().split("\t")
            # key: (chr, pos, ref, alt)
            key = (fields[0], fields[1], fields[2], fields[3])
            # value: (gene_name, pai3d_score, prediction)
            variant_dict[key] = (fields[4], float(fields[8]), fields[11])

    print(f"Dictionary built with {len(variant_dict)} entries")

    with open(pickle_path, "wb") as f:
        pickle.dump(variant_dict, f)

    print(f"Dictionary saved to pickle: {pickle_path}")

    return variant_dict


if __name__ == "__main__":
    # chromosome number passed as command line argument
    # e.g. python build_dict.py 22
    if len(sys.argv) != 2:
        print("Usage: python build_dict.py <chr_num>")
        sys.exit(1)

    chr_num = sys.argv[1]
    variant_dict = build_dict(chr_num)
    print(f"Ready. Example entry: {next(iter(variant_dict.items()))}")