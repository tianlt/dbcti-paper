
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dbcti-paper

<!-- badges: start -->

<!-- badges: end -->

This repository contains the code for the dbcti-paper. Working directory
should be set as dbcti-paper-master for running the codes.

## Data structure

  - `datasets` folder contains simulated data included in the paper.
    Real datasets are too large to be included and are available at the
    corresponding GEO datasets. The folder is also the directory used
    for store the generated data from the analysis codes.
  - `simulation_and_real_data` folder contains codes for generating
    datasets, and for DBCTI and other tools constructing cell
    trajectories. The file `dataset.R` contains code for subsetting the
    original GEO dataset (requried to be downloaded according to the GEO
    accession number listed in the data availability section in
    manuscript). File `simulated dataset.R` contains the code for
    generating the simulated dataset. For downstream analysis, the above
    two files should be run first. `plot_tools_function_upload.R`
    defines function of plotting result for manuscript. Other files
    contain codes for trajectory construction on datasets.
  - `evaluation` folder contains codes for evaluating performances of
    tools.
  - `case_study` folder contains codes for performing case study.

## File source

For files that are required as input in the code but not generated from
other codes, please refer to the data availability section in the
original paper.
