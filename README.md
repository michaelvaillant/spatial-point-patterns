# spatial-point-patterns

Reproducible spatial point-pattern analysis of **GEIPAN** UAP cases in Metropolitan France using non-stationary **Poisson point process models (PPM)**.  
The pipeline builds smoothed spatial covariates (demographic, environmental, infrastructure), produces **PDF figures** and **LaTeX tables**, and fits **spatstat** models with diagnostics.

## Contents

- `UAP_Analysis.R` — main analysis script (end-to-end)
- `data/` —  (see *Data & sources*)
- `cache/` — cached intermediate objects (created automatically)
- `output/figures/` — generated PDF figures
- `output/tables/` — generated LaTeX tables

## Requirements

- R (tested with R 4.x)
- Packages used by the script:
  - `sf`, `sp`, `spatstat`, `spdep`
  - `raster`, `png`, `pixmap`
  - `classInt`, `RColorBrewer`, `vcd`, `xtable`, `mnormt`

Install packages (example):

```r
install.packages(c(
  "sf","sp","spatstat","spdep","raster","png","pixmap",
  "classInt","RColorBrewer","vcd","xtable","mnormt"
))
