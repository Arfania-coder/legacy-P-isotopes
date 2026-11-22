# Changes

## 2.0.0

Version for the revised manuscript, "Oxygen isotopes in phosphate (δ18OP) reveal functional
separation between labile and stable legacy phosphorus in soil" (Biogeochemistry). It replaces
the scripts of v1.2.1, which accompanied the first submission.

- Data added, with a description of every column in `data/README.md`:
  `data/d18Op_hedley_fractions.csv` (48 observations), `data/precipitation_isotopes.csv`
  (22 samples) and `data/pullman_2nw_daily.csv` (daily station record, 2023–2024).
- One script per part of the paper (`R/01_data.R` to `R/10_seasonal_equilibrium.R`), run in
  order by `run_all.R`. They reproduce Tables 1 and 2, Supplementary Tables S1 to S3, the
  statistics of Sections 3.1 to 3.5 and of Supplementary Texts S1 and S2, and Figures 2 to 5,
  S1 and S2.
- Equilibrium δ18OP is calculated in the code (Chang and Blake, 2015) from the soil-water δ18O
  of each composite sample and the soil temperature measured on the sampling date. The
  principal component analysis and all figures use these values.
- Statistics follow Section 2.6 of the manuscript: additive linear models with type II analysis
  of variance, η², one-sample t-tests against equilibrium, and Tukey-adjusted comparisons of
  estimated marginal means.
- The soil-water line and the deuterium excess of soil water are calculated from the 12
  composite samples.
- The connectivity analysis is reproducible: the random-number generator and its seed are
  fixed, and the 4000 draws are written to `output/tables/Figure5_connectivity_draws.csv`.
- New: sensitivity of the equilibrium reference to the phosphate–water calibration
  (Supplementary Text S1, Table S3).
- New: `output/tables/manuscript_values_check.csv`, a line-by-line comparison of the values
  quoted in the manuscript with the script output.
- Figures are drawn at their final size (174 mm wide) without titles inside the figure, with
  one colour and one symbol per fraction in all of them.
- Title, authors and metadata updated (`CITATION.cff`, `.zenodo.json`); `LICENSE` added.
