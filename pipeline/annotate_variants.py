import gzip 
import csv 
import os
import sys

sys.path.append("/home/shatcher1/projects/PrimateAI-3D/pipeline")
from build_dict import build_dict

# paths

PAI3D_CHR22_PATH = "/home/shatcher1/projects/PrimateAI-3D/data/PrimateAI-3D.hg38.chr22.txt.gz"
PICKLE_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_variant_dict.pkl"
VCF_PATH = "/home/shatcher1/projects/PrimateAI-3D/data/1000g/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
OUTPUT_PATH = "/home/shatcher1/projects/PrimateAI-3D/pipeline/output/chr22_annotated.tsv"

# main function

def annotate_variants(vcf_path, output_path, variant_dict):
    
    with gzip.open(vcf_path, "rt") as vcf, open(output_path, "w", newline = "") as out:
        writer = csv.writer(out, delimiter = "\t")

        sample_ids = None 
        found = 0
        skipped = 0

        for line in vcf:
            
            # skip metadata lines -> any line that begins with ## in VCF format is metadata
            if line.startswith("##"):
                continue

            # header line - extract sample IDs and write output header
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t") # .strip() cleans up accidental spaces or hidden new lines
                                                  # .split("\t") looks across string for tabs -> "\t" , removes it, and cuts at that point
                                                  # EXAMPLE: #CHROM→POS→ID→REF→ALT→QUAL→FILTER→INFO→FORMAT→Sample1→Sample2 would turn into
                                                  # fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample1", "Sample2"]

                sample_ids = fields[9:]
                writer.writerow(["CHROM", "POS", "REF", "ALT", "PREDICTION"] + sample_ids)
                continue

            # data lines
            fields = line.strip().split("\t")
            key = (fields[0], fields[1], fields[3], fields[4])
            
            prediction = variant_dict.get(key)

            if prediction is None:
                skipped += 1
                continue

            genotypes = fields[9:]
            writer.writerow([fields[0], fields[1], fields[3], fields[4], prediction] + genotypes)
            found += 1
        
        print(f"Variants matched: {found}")
        print(f"Variants skipped: {skipped}")


if __name__ == "__main__":
    variant_dict = build_dict(PAI3D_CHR22_PATH, PICKLE_PATH)
    annotate_variants(VCF_PATH, OUTPUT_PATH, variant_dict)
    print(f"output written to {OUTPUT_PATH}")

