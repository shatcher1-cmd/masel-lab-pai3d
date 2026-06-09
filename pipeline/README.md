# PAI3D × 1000 Genomes Pipeline

## Overview
This pipeline connects PrimateAI-3D (PAI3D) pathogenicity scores with whole-genome 
sequencing data from the 1000 Genomes Project to test whether deleterious mutations 
are underdispersed across individuals — a prediction of Matheson & Masel 2023 arising 
from background selection at human-realistic mutation rates.

The current implementation runs on chromosome 22 as a proof of concept.

## Data Sources
**PrimateAI-3D (PAI3D)**
- Source: Illumina (Gao et al. 2023)
- File: `data/PrimateAI-3D.hg38.chr22.txt.gz`
- Contains pathogenicity predictions for missense variants genome-wide
- Variants scored on a continuous scale; binary prediction (benign/pathogenic) 
  used for current analysis
- Reference genome: GRCh38

**1000 Genomes Project**
- Source: Byrska-Bishop et al. 2021
- File: `data/1000g/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`
- 3,202 individuals across 26 populations, 30x coverage, phased
- Reference genome: GRCh38

## Pipeline Structure

build_dict.py → annotate_variants.py → compute_burden.py → compute_dispersion.py

### `build_dict.py`
Parses the PAI3D chr22 score file and builds a Python dictionary keyed by 
`(chromosome, position, ref, alt)` with PAI3D prediction as the value. 
Pickles the dictionary to disk for reuse across runs.

- Input: `data/PrimateAI-3D.hg38.chr22.txt.gz`
- Output: `pipeline/output/chr22_variant_dict.pkl`
- Entries: 1,476,838 scored variants on chr22

### `annotate_variants.py`
For each variant in the 1000G VCF, looks up its PAI3D prediction in the dictionary.
Variants with no PAI3D score (non-coding variants, INDELs) are skipped. Includes 
complement strand lookup to handle strand orientation differences between datasets.

- Input: 1000G chr22 VCF + PAI3D dictionary
- Output: `pipeline/output/chr22_annotated.tsv`
- Columns: CHROM, POS, REF, ALT, PREDICTION, [3202 individual genotype columns]
- Matched variants: 7,478 (633 pathogenic, 6,845 benign)
- Skipped variants: 1,059,079 (non-coding, INDELs, unscored)

### `compute_burden.py`
For each individual, counts the total number of pathogenic and benign variant 
copies they carry. Counts actual allele copies (0, 1, or 2) rather than binary 
presence/absence, reflecting the biological difference between heterozygous and 
homozygous carriers.

- Input: `pipeline/output/chr22_annotated.tsv`
- Output: `pipeline/output/chr22_burden.tsv`
- Columns: SAMPLE_ID, PATHOGENIC_BURDEN, BENIGN_BURDEN

### `compute_dispersion.py`
Computes the index of dispersion (variance / mean) for pathogenic and benign 
burden counts across all 3,202 individuals. Underdispersion (index < 1) indicates 
individuals vary less in mutation burden than expected under independence.

- Input: `pipeline/output/chr22_burden.tsv`
- Output: `pipeline/output/chr22_dispersion.tsv`

## Preliminary Results (chr22)

| Category | Mean | Variance | Index of Dispersion |
|---|---|---|---|
| Pathogenic | 6.17 | 5.29 | **0.857** |
| Benign | 300.09 | 2405.07 | **8.014** |

Pathogenic variants are underdispersed (IoD < 1), consistent with the theoretical 
prediction of Matheson & Masel 2023. Benign variants are overdispersed as expected 
given the absence of purifying selection and the influence of population structure.

## Environment
- Conda environment: `pai3d_public`
- Reference genome: GRCh38 throughout
- Currently running locally; HPC scaling to be determined

## Known Limitations
- chr22 only — full genome scaling in progress
- Binary PAI3D prediction used; continuous scores to be incorporated
- Population stratification not yet implemented
- Pipeline output to be refactored into a structured data object

## Next Steps
- Refactor annotated output into a structured data object including raw PAI3D scores
  and additional variant metadata
- Scale pipeline to full genome
- Add population stratification analysis
- Implement LD curve analysis (Good 2022 framework)
