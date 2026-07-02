# CLAUDE.md — Project Context for AI Assistance

## Project Overview
This pipeline connects PrimateAI-3D (PAI3D) pathogenicity scores with 1000 Genomes 
Project whole-genome sequencing data to test whether deleterious mutations are 
underdispersed across individuals — a prediction of Matheson & Masel 2023 arising 
from background selection. This is undergraduate research in the Masel Lab 
(PI: Joanna Masel, Department of Ecology and Evolutionary Biology, University of Arizona).

## Environment
- Conda environment: pai3d_public
- OS: WSL2 Ubuntu on Windows
- Python 3.x, numpy 1.23.5
- VS Code with WSL2
- Reference genome: GRCh38 throughout
- Local laptop for now; HPC scaling TBD

## Repository Structure
- data/PrimateAI-3D.hg38.chr22.txt.gz — PAI3D scores, chr22, GRCh38
- data/PrimateAI-3D.hg38.chr22.txt.gz.tbi — tabix index
- data/1000g/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
- data/1000g/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
- data/1000g/20130606_g1k_3202_samples_ped_population.txt — population metadata (see below)
- pipeline/README.md
- pipeline/build_dict.py
- pipeline/variants.py          — Variants class (renamed from variant_record.py / VariantRecord)
- pipeline/individual.py        — Individual class
- pipeline/annotate_variants.py
- pipeline/output/chr22_variant_dict.pkl
- pipeline/output/chr22_individuals.pkl  — pickled list of Individual objects (primary output)

## PAI3D Data Structure
File: PrimateAI-3D.hg38.chr22.txt.gz
Columns (0-indexed):
- 0: chr
- 1: pos
- 2: non_flipped_ref
- 3: non_flipped_alt
- 4: gene_name (NOTE: appears to store transcript IDs e.g. ENST00000252835.5, not gene symbols — still unclarified with Joanna)
- 5: change_position_1based
- 6: ref_aa
- 7: alt_aa
- 8: score_PAI3D (continuous, 0-1)
- 9: percentile_PAI3D
- 10: refseq
- 11: prediction (benign/pathogenic, PAI3D's own threshold = 0.821; now redundant since we derive prediction ourselves — candidate for removal from dict value)

Dictionary key: (fields[0], fields[1], fields[2], fields[3])
Dictionary value: (fields[4], float(fields[8]), fields[11]) — (gene_name, pai3d_score, prediction)
1,476,838 entries for chr22

## 1000G Data Structure

### VCF
File: 1kGP chr22 VCF, Byrska-Bishop et al. 2021, 30x coverage, GRCh38, 3,202 individuals
Source: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
VCF columns (0-indexed):
- 0: CHROM (chr22)
- 1: POS
- 2: ID
- 3: REF
- 4: ALT
- 5: QUAL
- 6: FILTER
- 7: INFO
- 8: FORMAT (GT = genotype)
- 9+: individual genotypes (0|0, 0|1, 1|0, 1|1)

### Population metadata (ped file)
File: 20130606_g1k_3202_samples_ped_population.txt
Source: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
Format: space-delimited, one row per sample, NOT tab-delimited like the VCF
Columns: FamilyID, SampleID, FatherID, MotherID, Sex, Population, Superpopulation
Only SampleID (col 1) and Population (col 5) are currently used; FatherID/MotherID 
(pedigree info) are not used by this pipeline
All 3,202 samples in the VCF matched successfully to this file — verified, 
zero individuals with population=None after integration

## Class Structure

### Variants (variants.py)
Renamed from VariantRecord / variant_record.py.
Represents a single variant for a single individual. Only created when genotype_copies > 0.
Fields:
- chromosome: str — e.g. "chr22"
- gene: str — from PAI3D fields[4] (transcript ID, not gene symbol)
- position: int — genomic position
- reference_allele: str
- alternate_allele: str
- genotype_copies: int — 0, 1, or 2 alternate allele copies
- pai3d_score: float — raw continuous score from PAI3D fields[8]
- pai3d_prediction: str — "pathogenic" if pai3d_score > 0.821 else "benign" (computed, not stored from file)

### Individual (individual.py)
Represents one person from the 1000G dataset.
Fields:
- sample_id: str — e.g. "HG00096"
- population: str — now populated from ped file (e.g. "GBR"); no longer defaults to None 
  in practice, though the parameter default remains None for flexibility
- variants: list of Variants objects

Methods:
- add_variant(variant) — appends a Variants object to variants list
- compute_burden(threshold=0.821) — returns dict {"pathogenic": int, "benign": int}
- compute_dispersion(individuals, threshold=0.821) — takes list of Individual objects, 
  returns dict {"pathogenic": float, "benign": float}
- compute_dispersion_by_population(individuals, population, threshold=0.821) — filters 
  individuals list to matching population label first, then calls compute_dispersion; 
  returns None if no individuals match the given population

## Pipeline Flow
1. build_dict.py — parses PAI3D file, pickled dict keyed by (chr, pos, ref, alt),
   value is (gene_name, pai3d_score, prediction) tuple
2. annotate_variants.py — loads population map from ped file, joins 1000G VCF with 
   PAI3D dict, builds Individual objects with population assigned, pickles list of 
   3,202 Individual objects to chr22_individuals.pkl
3. compute_burden(), compute_dispersion(), and compute_dispersion_by_population() are 
   methods on Individual, not standalone scripts

## Preliminary Results (chr22)
- Matched variants: 7,478
- Skipped: 1,059,079 (non-coding, INDELs, unscored — expected)
- Example individual HG00096: 193 variants, pathogenic burden 2, benign burden 278
- Cohort-wide dispersion (class-based pipeline, all 3,202 individuals):
  - Pathogenic: 0.9576
  - Benign: 8.1298
- NOTE: pathogenic dispersion (0.9576) differs from original TSV-based pipeline result 
  (0.857) — benign dispersion is consistent (8.13 vs 8.014, negligible). Suspected cause: 
  threshold derivation method (own 0.821 cutoff vs. PAI3D's internal fields[11] label) — 
  not yet root-caused
- Population-stratified dispersion not yet computed/reported (method exists, untested 
  with real population data as of this writing)

## Key Decisions Made
- Per-chromosome dictionary organization
- Raw PAI3D score stored (fields[8]); prediction derived at threshold 0.821
- Genotype counting: actual copies (0, 1, or 2) not binary presence/absence
- Variants only created when genotype_copies > 0
- One Individual object per sample; variants held as list of Variants objects
- Primary output: pickled list of Individual objects (chr22_individuals.pkl)
- Tabular export to be added later via to_dataframe() method
- Complement strand lookup included in annotate_variants.py
- Standalone scripts for now; Individual and Variants imported as modules
- VariantRecord renamed to Variants (and file renamed variant_record.py -> variants.py) 
  per Joanna's request
- Population now integrated via ped file; population stratification supported via 
  compute_dispersion_by_population()

## Open Questions
1. fields[4] in PAI3D stores transcript IDs (ENST...) not gene symbols — is this 
   the correct field for gene name, or should a different field be used?
2. Why does pathogenic dispersion differ between old pipeline (0.857) and new pipeline 
   (0.9576) when benign dispersion is consistent? Investigate threshold/prediction logic.
3. fields[11] (PAI3D's own prediction label) is redundant in the dictionary since 
   prediction is now derived independently — candidate for removal from build_dict.py
4. Next analysis priorities: population stratification (now possible), LD curves, 
   gene-level stratification?
5. Should this pipeline move to its own repository?
6. Multi-chromosome scaling: code refactor is straightforward (loop instead of hardcode), 
   but full 1000G download size and storage likely requires HPC rather than laptop — 
   raise with Joanna before downloading remaining chromosomes

## Future Datasets
- Quebec Biobank — potential future replacement or complement to 1000G data
  (access requirements and data format TBD)

## Coding Preferences
- snake_case naming
- Comments explaining what each block does
- Understand reasoning before implementation
- Flag low confidence rather than guess
- Imports at the top of every file, never mid-script
- Readability preferred over efficiency (explicit loops over list comprehensions 
  where it aids clarity)
- Standalone scripts for now; modular refactor later