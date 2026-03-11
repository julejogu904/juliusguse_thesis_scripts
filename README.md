# Bachelor's Thesis Repository

This repository contains the scripts, environment specifications and intermediate, output and final files from the phylogenetics pipeline.

## Contents

The repository includes:

- custom scripts developed and used for the pipeline
- the `environment.yml` file for reproducing the software setup
- relevant intermediate files (s4predout, parsed_s4pred, SR3_recoded, indel_transfer, trimmed)
- concatenated alignments for all final datasets
- final phylogenetic tree files

## Data types

The analyses were based on four final datasets:

- amino acid sequence data (einsi)
- SR3-recoded data
- S4PRED predicted secondary structure data with threshold 1
- S4PRED predicted secondary structure data with threshold 5

## Reproducibility
 
All files required to reproduce the main analyses of this thesis are provided in this repository.

## Notes

Some external software tools used in the workflow include WhereDoGGo, IQ-TREE, iTOL, PhyKIT and S4PRED.  
Please create the conda environment from `environment.yml` before running the scripts.
