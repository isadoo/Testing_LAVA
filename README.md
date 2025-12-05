# Testing LAVA method:
# Simulation-Based Method Comparison for Detecting Selection in Quantitative Traits

This repository contains simulation scripts, data preparation and analysis workflows, and results from comparing different methods of detecting selection on quantitative traits in structured populations.

## Project Overview

This study evaluates three different statistical methods for detecting selection on quantitative traits using simulated data. We simulated populations with varying structures and selection regimes using Quantinemo, then compared the performance of three methods (including a variation of LAVA which considers habitat information):

1. **Qst-Fst** approach
2. **Driftsel** (From Ovaskainen and collaborators 2011 & 2013)
3. **LAVA** (Log-Ratio of Ancestral Variances - from do O and collaborators 2025)
4. **LAVA w/ habitat**

## Workflow

### 1. Simulations (Quantinemo)
We ran individual-based simulations using Quantinemo with:
- **Population structures**: 
  - Stepping Stone (SS) - 18 populations
  - Island Model (IM) - 18 and 9 populations
  - Hierarchical - 18 populations
- **Selection regimes**: Varying selection coefficients (wdiff), environmental variance (wvar), and correlation between selection and relatedness (correlation).

### 2. Data Preparation
All simulated data were standardized for analysis across methods:
- Generated F1 individuals through controlled breeding designs
- Calculated phenotypes for F1 generation
- Estimated coancestry matrices
- Ensured identical datasets for all three methods

### 3. Method Application
Applied each of the three methods to the exact same prepared datasets to obtain:
- P-values (Qst-Fst and Driftsel)
- S-values (LAVA)

### 4. Performance Evaluation
Built ROC curves and calculated AUC (Area Under the Curve) to compare method performance across different scenarios.

## Repository Structure

```
Testing_LAVA/
├── simulation_scripts/          # Quantinemo configuration files (.ini) for running simulations
│                                 # Includes all population structures (SS, IM, Hierarchical) with
│                                 # varying selection coefficients and environmental variance
│
├── method_testing_scripts/      # Scripts for data preparation and method application
│   ├── breeding_tests/          # R scripts for creating F1 generations under different
│   │                            # breeding designs (varying pop numbers, F1 numbers, sire/dam ratios)
│   └── testing_scripts/         # R scripts implementing each method (Qst-Fst, LAVA, Driftsel)
│                                # and combining results across simulations
│
├── tools/                       # Utility R functions used across all analyses
│                                # (coancestry calculations, F1 creation, method implementations)
│
├── raw_data/                    # Combined results (p-values/s-values) from all method applications
│                                # TSV files containing test statistics for each locus and scenario
│
├── analysis_scripts/            # R scripts for generating ROC curves and calculating performance
│                                # metrics (AUC, TPR at specific FPR thresholds)
│
└── results/                     # Final results: AUC summaries (CSV), ROC 
                                 # and performance tables (CSV/TEX)
```

### Directory Descriptions

**simulation_scripts/**  
Contains all Quantinemo `.ini` configuration files used to generate simulated datasets. These define population parameters (size, structure, migration rates), selection coefficients (wdiff), environmental variance (wvar), and genetic architecture for quantitative traits.

**method_testing_scripts/**  
Two subdirectories:
- `breeding_tests/`: Scripts that create F1 individuals from simulated parental populations using various breeding designs - `testing_data_preparation.r` is the script we used on the basic breeding desing, which we at times refer to as "full breeding desing".
- `testing_scripts/`: Scripts that apply Qst-Fst, LAVA, LAVA w/ habitat, and Driftsel methods to the prepared data, ensuring all methods analyze identical datasets

**tools/**  
Reusable R functions including:
- Coancestry matrix estimation from genotype data
- F1 generation from parental populations
- Driftsel implementation
- LAVA implementation (with and without habitat information)
- Kinship calculations from pedigree data (unused in final versions)

**raw_data/**  
Tab-separated files containing the raw test statistics (p-values for Qst-Fst and LAVA, s-values for Driftsel) for each trait across all simulated scenarios and replicates.

**analysis_scripts/**  
R scripts that generate ROC curves by comparing TPR vs FPR. Calculates Area Under the Curve (AUC) and True Positive Rates at specified False Positive Rate thresholds.

**results/**  
Final outputs including:
- AUC summary tables comparing method performance across scenarios
- ROC curve plots
- Performance metrics tables in CSV and LaTeX formats
