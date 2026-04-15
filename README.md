README: NCA in Psychological Data
================

This repository contains code and materials for a simulation study on
**Necessary Condition Analysis (NCA)** in psychological research.

## Goal

Assess how much NCA estimates change when the latent relationship is
held constant, but typical psychometric features of observed data
vary.

## Key Idea

NCA operates on the **observed score distribution**, not on latent
constructs.

In psychological data, observed scores depend on: 

- measurement error (reliability) 
- ordinal scaling (Likert categories) 
- threshold placement
- skewness
- sample selection
- ...

The project tests how strongly these features affect NCA results.

## Design

Paired simulation:

1.  Generate latent variables (fixed relation)
2.  Compute **latent-reference NCA**
3.  Create observed datasets from the same latent data
4.  Recompute NCA under different conditions
5.  Compare estimates

## Simulation Factors

- **Reliability** (factor loadings)
- **Ordinal scaling** (thresholds, category structure)
- **Score construction** (sum scores)
- **Range restriction**

## Main Files

- `simulations/nca_unified_paired_single_script.R` — main simulation script  
- `simulations/simulation_results.csv` — raw results  
- `simulations/results_table.csv` — aggregated results  
- `simulations/simulation_results_n.csv` — additional sample size results 
- `simulations/results_table_n.csv` — additional sample size aggregated results  
- `nca_paper.qmd` — manuscript
- `nca_supplementary.cmq` - supplementary materials

## Reproducibility

Uncomment and run:

``` r
source("nca_unified_paired_single_script.R")
```
