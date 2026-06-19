# Necessary Condition Analysis project

This repository contains the data, simulation code, manuscript sources, and rendered outputs for the NCA psychometrics project.

## Repository map

- `data/`: Scopus export data used in the manuscript figures.
  - `nca_all.csv`: publications matching the broader NCA search.
  - `nca_psic.csv`: psychology subset.
- `simulations/`: R scripts and simulation outputs.
  - `nca_unified_paired_single_script_rev.R`: full paired simulation framework.
  - `nca_simulation_fast.R`: faster CE-FDH implementation and main simulation run.
  - `nca_script_concordance.R`: concordance check between the fast implementation and `NCA`.
  - `simulation_results_*.csv`: replication-level outputs.
  - `results_table_*.csv`: aggregated summaries used by the paper/supplement.
- `R/`: plotting helpers and exploratory plotting code.
  - `plots.R`: ggplot summaries for simulation result tables.
- `paper/`: manuscript, supplement, references, figures, and rendered outputs.
  - `nca_paper.qmd`: main Quarto manuscript source.
  - `nca_supplementary.qmd`: supplementary materials source.
  - `nca_refs.bib`: bibliography.
  - `nca_paper.pdf`, `nca_paper.html`, `nca_paper.docx`: rendered manuscript outputs.
  - `nca_supplementary.html`: rendered supplement.
  - `_extensions/wjschne/apaquarto/`: local Quarto APA extension files.

## Main outputs

- Main simulation summary: `simulations/results_table_fast_10k.csv`
- Main simulation replications: `simulations/simulation_results_fast_10k.csv`
- Sample-size summary: `simulations/results_table_n.csv`
- Sample-size replications: `simulations/simulation_results_n.csv`
- Fast/NCA concordance output: `simulations/simulation_results_concordance.csv`

Large generated files are kept in the repository because the manuscript and supplement read from the saved result tables.

## Reproducing

Open `necessary-condition-analysis.Rproj` in RStudio, or use the repository root as the working directory.

Install the R packages used across the project:

```r
install.packages(c(
  "dplyr", "purrr", "tibble", "tidyr", "ggplot2",
  "lavaan", "NCA", "furrr", "future", "patchwork",
  "here", "knitr"
))
```

Render the manuscript and supplement with Quarto:

```sh
quarto render paper/nca_paper.qmd
quarto render paper/nca_supplementary.qmd
```

Simulation scripts include run blocks and parallel worker settings near the bottom of each file. Check those parameters before sourcing them, especially on local machines.

## License

See `LICENSE`.
