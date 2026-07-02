import gzip 
import pickle
import os

# paths

PAI3D_CHR22_PATH = "/home/shatcher1/projects/PrimateAI-3D/data/PrimateAI-3D.hg38.chr22.txt.gz"
PICKLE_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_variant_dict.pkl"

# Load and return the cached dictionary if the pickle file 
# already exists

def build_dict(pai3d_path, pickle_path):
    if os.path.exists(pickle_path):
        print("Loading existing dictionary from pickle...")
        with open(pickle_path, "rb") as f:
            return pickle.load(f)

    print("Building dictionary from PAI3D file...") 
    variant_dict = {}

    with gzip.open(pai3d_path, "rt") as f:
        next(f) # skip's header
        for line in f:
            fields = line.strip().split("\t")
            key = (fields[0], fields[1], fields[2], fields[3]) # uses 'chr' , 'pos' , 'non_flipped_ref' , 'non_flipped_alt' as the key -> fields[0:3]
            variant_dict[key] = (fields[4], float(fields[8]), fields[11])   # value is a tuple of (gene_name, pai3d_score, prediction) — fields[4], [8], [11]


    print(f"Dictionary built with {len(variant_dict)} entries")

    with open(pickle_path, "wb") as f:
        pickle.dump(variant_dict, f)
    print("Dictionary saved to pickle")

    return variant_dict
    

if __name__ == "__main__":
    variant_dict = build_dict(PAI3D_CHR22_PATH, PICKLE_PATH)
    print(f"Ready. Example entry: {next(iter(variant_dict.items()))}")

