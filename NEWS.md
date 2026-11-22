# Changes

## 2.1.0 (2026-11-22)

Completes the archive for the revised manuscript. Version 2.0.0 held the concentration and
isotope statistics and Figure 3; this version adds the remaining analyses and their data, so
that `run_all.R` reproduces every table, every quoted statistic and every data figure of the
paper.

- Data: `data/d18Op_hedley_fractions.csv` now has the precision of the laboratory data sheet
  and the soil-water δ2H. `data/precipitation_isotopes.csv` (22 samples) and
  `data/pullman_2nw_daily.csv` (daily station record, 2023–2024) are new. The daily record
  replaces `data/soil_temperature_summary.csv`; Supplementary Table S2 is calculated from it.
- New scripts: `02_station_climate.R` (Table S2, Figure S2), `04_water_isotopes.R` (Figure 2),
  `07_pca.R` (Figure 4), `08_connectivity.R` (Figure 5) and `10_seasonal_equilibrium.R`
  (Text S2, Figure S1). The scripts of 2.0.0 are renumbered `03`, `05`, `06` and `09`.
- The principal component analysis and all figures use equilibrium values calculated with
  Eq. 1 (Chang and Blake, 2015).
- The soil-water line and the deuterium excess of soil water are calculated from the 12
  composite samples.
- The connectivity analysis is reproducible: the random-number generator and its seed are
  fixed, and the 4000 draws are written to `output/tables/Figure5_connectivity_draws.csv`.
- Figures are drawn at their final size (174 mm wide) without titles inside the figure, with
  one colour and one symbol per fraction in all of them.
- `output/tables/manuscript_values_check.csv` now covers 263 quoted values.

## 2.0.0 (2026-11-22)

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
