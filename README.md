# legacy-P-isotopes

Data and R code for

> Arfania H, Kayler Z, Strawn D, Brooks E, Feber G, Laan M. Oxygen isotopes in phosphate
> (δ<sup>18</sup>O<sub>P</sub>) reveal functional separation between labile and stable legacy
> phosphorus in soil. *Biogeochemistry* (in revision).

The repository holds the measurements behind the paper and the R scripts that turn them into
its tables, statistics and Figure 3. The measurements are the oxygen isotope composition of
phosphate (δ<sup>18</sup>O<sub>P</sub>) and the P concentration of four sequentially extracted
(Hedley) fractions in 12 composite soil samples from the R. J. Cook Agronomy Farm (Pullman,
Washington, USA), with the soil-water δ<sup>18</sup>O and soil temperature needed to calculate
isotopic equilibrium.

## Contents

| Path | Content |
| --- | --- |
| `data/d18Op_hedley_fractions.csv` | 48 observations: 4 fractions × 2 landscape positions × 2 tillage systems × 3 sampling dates |
| `data/soil_temperature_summary.csv` | Soil temperature at the Pullman 2 NW station by period and depth (Supplementary Table S2) |
| `data/README.md` | Description of every column, with units and sources |
| `R/00_setup.R` | Constants, equilibrium equations, statistics and formatting helpers |
| `R/01_data.R` | Reads the data and calculates equilibrium δ<sup>18</sup>O<sub>P</sub> and the departure from it |
| `R/02_fractionation.R` | P concentrations of the fractions: Table 1 and Section 3.1 |
| `R/03_isotopes.R` | δ<sup>18</sup>O<sub>P</sub> of the fractions: Table 2 and Sections 3.3 to 3.5 |
| `R/04_equilibrium_sensitivity.R` | Sensitivity of the equilibrium reference: Supplementary Text S1 and Table S3 |
| `R/05_figure3.R` | Figure 3 |
| `run_all.R` | Runs the scripts in order and compares the results with the values quoted in the paper |
| `output/` | Tables and figure written by the scripts |

## What the scripts reproduce

| In the paper | Script | Output in `output/` |
| --- | --- | --- |
| Supplementary Table S1 (equilibrium and Δδ<sup>18</sup>O<sub>P</sub> of the 48 observations) | `01_data.R` | `tables/TableS1_d18Op_equilibrium.csv` |
| Table 1 and Section 3.1 (factor means, ANOVA, position × tillage means, P stocks) | `02_fractionation.R` | `tables/Table1_*.csv`, `tables/Sect3.1_*.csv` |
| Table 2 and Sections 3.3 and 3.4 (means, t-tests, ANOVA, η², interactions, Tukey comparisons) | `03_isotopes.R` | `tables/Table2_*.csv`, `tables/Sect3.4_*.csv` |
| Pairwise isotopic separation between fractions (Section 3.5, Figure 5c) | `03_isotopes.R` | `tables/Figure5c_pairwise_separation.csv` |
| Supplementary Text S1 and Table S3 (four phosphate–water calibrations, extrapolation below 3 °C) | `04_equilibrium_sensitivity.R` | `tables/TableS3_calibration_sensitivity.csv`, `tables/TextS1_extrapolation_and_range.csv` |
| Figure 3 | `05_figure3.R` | `figures/Figure3_d18Op_vs_equilibrium.png` and `.pdf`, `tables/Figure3_equilibrium_bands.csv` |

Not part of this version: the principal component analysis (Figure 4), the Bayesian
isotopic-connectivity model (Figure 5a, b and d), the water-isotope figure (Figure 2) and the
supplementary figures (Figures S1 and S2). They will be added together with the data they
need: the soil-water δ<sup>2</sup>H values, the precipitation isotope data and the daily
station record.

## Running the code

Requirements: R 4.1 or later with the packages `emmeans` and `ggplot2` (3.4.0 or later). The R
and package versions used for the archived outputs are recorded in `output/sessionInfo.txt`.

```r
install.packages(c("emmeans", "ggplot2"))
```

Run everything from the repository root, either in a terminal

```
Rscript run_all.R
```

or in R or RStudio (open `legacy-P-isotopes.Rproj`, or set the working directory to this folder)

```r
source("run_all.R")
```

The run takes a few seconds. It writes the tables to `output/tables/`, Figure 3 to
`output/figures/` and the R session information to `output/sessionInfo.txt`.

## Checking the paper against the code

Each statistic quoted in the paper that these scripts cover is registered in the script that
calculates it, together with the value printed in the paper. `run_all.R` collects them (175
values) in `output/tables/manuscript_values_check.csv`, with one row per value: where it is
quoted, the value in the paper, the value calculated here formatted to the same number of
decimals, whether the two agree, and a note where they do not. The run ends with a summary of
this comparison.

## Methods in brief

- **Equilibrium.** Equilibrium δ<sup>18</sup>O<sub>P</sub> is calculated with the
  phosphate–water equation of Chang and Blake (2015), 1000 ln α = 14.43 × 10<sup>3</sup>/T − 26.54
  (T in kelvin), from the soil-water δ<sup>18</sup>O of each composite sample and the soil
  temperature at 10 cm measured on the sampling date. Δδ<sup>18</sup>O<sub>P</sub> is the measured
  value minus the equilibrium value.
- **Linear models.** δ<sup>18</sup>O<sub>P</sub> and Δδ<sup>18</sup>O<sub>P</sub> are modelled
  with fraction, landscape position, tillage and sampling date as additive fixed factors
  (type II analysis of variance; η² is the factor sum of squares divided by the total sum of
  squares). Interactions with fraction are tested by adding one term at a time. Concentrations
  are modelled with position, tillage and date, and the position × tillage interaction is
  tested by adding it to that model.
- **Comparisons.** Each fraction is tested against equilibrium with a one-sample t-test of
  Δδ<sup>18</sup>O<sub>P</sub> against zero (n = 12). Fractions and sampling dates are compared
  with Tukey-adjusted contrasts of estimated marginal means (`emmeans`).
- **Equilibrium bands of Figure 3.** For each composite sample, equilibrium is evaluated over
  the range of daily mean soil temperature at 5 to 20 cm in its sampling month.
- **Sensitivity.** All 48 observations are recalculated with the calibrations of Longinelli and
  Nuti (1973), Lécuyer et al. (2013) and Pucéat et al. (2010).

## Data

The data files are described in [`data/README.md`](data/README.md). The δ<sup>18</sup>O<sub>P</sub>,
soil-water δ<sup>18</sup>O and soil-temperature values are those of Supplementary Table S1, the
concentrations those of Table 1, and the station summary that of Supplementary Table S2.

## Citation

Please cite the archived version of this repository on Zenodo
(<https://doi.org/10.5281/zenodo.18752516>, which resolves to the latest version) and the
article above. Citation metadata are in `CITATION.cff`.

## License

MIT License; see `LICENSE`.
