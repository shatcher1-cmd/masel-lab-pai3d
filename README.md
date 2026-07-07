# masel-lab-pai3d
A bioinformatics pipeline built for undergraduate research in the Masel Lab 
(Department of Ecology and Evolutionary Biology, University of Arizona, PI: Joanna Masel).

## What this project does
This pipeline connects two large genomic datasets to test a prediction from 
evolutionary biology: that deleterious (harmful) mutations should be more evenly 
distributed across individuals than you would expect by chance — a pattern called 
underdispersion.

To test this, the pipeline:
1. Loads pathogenicity scores from PrimateAI-3D, a machine learning model trained 
   on primate genomes that scores how harmful each genetic variant is likely to be
2. Matches those scores against whole-genome sequencing data from 3,202 individuals 
   in the 1000 Genomes Project
3. Builds a per-individual data object storing every variant that person carries 
   along with its pathogenicity score
4. Computes the index of dispersion (variance / mean) across the cohort to measure 
   how evenly mutations are distributed

## Data sources
- **PrimateAI-3D** (Gao et al. 2023) — variant pathogenicity scores, GRCh38
- **1000 Genomes Project** high-coverage WGS (Byrska-Bishop et al. 2021) — 3,202 
  individuals, 30x coverage, GRCh38
- **1000G population metadata** — sample-level ancestry labels (26 population codes 
  across 5 superpopulations)

## Pipeline structure
```
pipeline/
    build_dict.py          — parses PAI3D scores into a fast lookup dictionary (accepts chr_num argument)
    annotate_variants.py   — matches 1000G variants against PAI3D scores, builds per-individual 
                             data objects, assigns population labels (accepts chr_num argument)
    individual.py          — Individual class: holds each person's variant list and methods 
                             for burden and dispersion analysis
    variant.py             — Variant class: stores per-variant data fields
    iod_analysis.ipynb     — analysis notebook: IoD table and plots vs threshold, chromosome, population
```

## Scaling
The pipeline runs on UA's Puma HPC cluster via a SLURM job array that processes 
all 22 autosomes in parallel. Each chromosome job:
1. Builds a PAI3D lookup dictionary from the per-chromosome score file
2. Annotates all 3,202 individuals with their variants for that chromosome
3. Pickles the resulting Individual objects to disk

## Preliminary results
- Chromosomes processed: 1-22 (autosomes only; chrX excluded pending sex-aware genotype handling)
- Individuals: 3,202 per chromosome
- Chr22 pathogenic index of dispersion: ~0.96 (below 1, consistent with underdispersion)
- Chr22 benign index of dispersion: ~8.13 (overdispersed, as expected)
- Per-chromosome pathogenic IoD ranges approximately 0.65-1.22 at the default threshold (0.821)

These results are preliminary and analysis is ongoing.

## Environment
- Python 3.x, numpy, pandas, matplotlib
- Conda environment: pai3d_public
- WSL2 Ubuntu (local dev) / UA Puma HPC cluster (scaling)
- Reference genome: GRCh38 throughout

## Status
All 22 autosomes complete. Analysis notebook in development — IoD as a function of 
threshold, chromosome, and population.
