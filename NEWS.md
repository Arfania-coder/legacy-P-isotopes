# Changes

## 2.0.0

Version for the revised manuscript, "Oxygen isotopes in phosphate (δ18OP) reveal functional
separation between labile and stable legacy phosphorus in soil" (Biogeochemistry). It replaces
the scripts of v1.2.1, which accompanied the first submission.

- Data added: `data/d18Op_hedley_fractions.csv` (48 observations) and
  `data/soil_temperature_summary.csv`, with a description of every column in `data/README.md`.
- Equilibrium δ18OP is calculated in the code (Chang and Blake, 2015) from the soil-water δ18O
  of each composite sample and the soil temperature measured on the sampling date.
- Statistics follow Section 2.6 of the manuscript: additive linear models with type II analysis
  of variance, η², one-sample t-tests against equilibrium, and Tukey-adjusted comparisons of
  estimated marginal means.
- New: sensitivity of the equilibrium reference to the phosphate–water calibration
  (Supplementary Text S1, Table S3).
- New: `output/tables/manuscript_values_check.csv`, a line-by-line comparison of the values
  quoted in the manuscript with the script output.
- Figure 3 is drawn at its final size (174 mm) without a title inside the figure.
- Title, authors and metadata updated (`CITATION.cff`, `.zenodo.json`); `LICENSE` added.
