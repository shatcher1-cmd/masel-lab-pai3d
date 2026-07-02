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
3. Builds a per-individual data object storing every harmful variant that person carries
4. Computes the index of dispersion (variance / mean) across the cohort to measure 
   how evenly mutations are distributed

## Data sources

- **PrimateAI-3D** (Gao et al. 2023) — variant pathogenicity scores, GRCh38
- **1000 Genomes Project** high-coverage WGS (Byrska-Bishop et al. 2021) — 3,202 
  individuals, 30x coverage, GRCh38
- **1000G population metadata** — sample-level ancestry labels (26 population codes 
  across 5 superpopulations)

## Pipeline structure
build_dict.py          → parses PAI3D scores into a fast lookup dictionary
annotate_variants.py   → matches 1000G variants against PAI3D scores, builds per-individual data objects, assigns population labels
individual.py          → Individual class: holds each person's variant list and methods for burden and dispersion analysis
variants.py            → Variants class: stores per-variant data fields


## Preliminary results (chr22)

- 7,478 variants matched to PAI3D scores across 3,202 individuals
- Pathogenic index of dispersion: **0.96** (below 1, consistent with underdispersion)
- Benign index of dispersion: **8.13** (overdispersed, as expected)

These results are preliminary. The 1000 Genomes Project dataset is being used 
as a starting point for pipeline development and the dataset will be updated in 
a future iteration of this project.

## Status

Chr22 complete. Currently scaling to all chromosomes via SLURM on UA's Puma 
HPC cluster.

## Environment

- Python 3.x, numpy, conda
- WSL2 Ubuntu / UA Puma HPC cluster
- Reference genome: GRCh38 throughout
