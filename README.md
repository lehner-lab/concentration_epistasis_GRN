# concentration_epistasis_GRN
# Overview

Welcome to the GitHub repository for concentration_epistasis_GRN.

# Brief description of the study: 
Both the effects of mutations and how they interact can change across conditions. Mutations can have diverse biochemical effects, for example having strong or weak effects on the stability of proteins, but also changing their expression and binding to other molecules. Here, we use a well-established thermodynamic model of gene regulation by one of the best-understood proteins, the CI transcriptional repressor of phage lambda, to better understand how diverse mutations interact and how both mutational effects and the interactions between mutations change in response to a change in the expression level of the protein itself. 

# Required Software

To run the R codes, the following software and associated packages are needed:

* **[R](https://www.r-project.org/) >=v3.5.2** (ggplot2, ggpubr, reshape2, rootSolve, seqinr, stringr)

# Codes
## file 1: Function_thermodynamics_model.R

Functions to calculate phenotypic outputs based on mutations affecting each parameter or combination of parameters, based on thermodynamics model of CI protein folding and CI-OR binding. 

## file 2: making_data_frames.R

This file contains codes to generate data modelling mutational effects. 

## file 3: plotting.R

Plotting, using the data lists generated with file 2 (making_data_frames.R). 

