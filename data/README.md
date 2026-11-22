# Data

Three files. Missing values are empty fields. Isotope ratios are in ‰ relative to VSMOW.

## `d18Op_hedley_fractions.csv`

One row per fraction and composite soil sample (48 rows). On each of three winter dates, one
composite sample (12 subsamples, 0–15 cm) was taken from each of four plots (two landscape
positions in each of two tillage fields) at the R. J. Cook Agronomy Farm near Pullman,
Washington, USA. Each of the 12 composite samples was extracted sequentially into four P
fractions.

| Column | Unit | Description |
| --- | --- | --- |
| `tillage` | – | `CT` = conventional tillage, `NT` = no-till |
| `position` | – | `Top` = top slope (erosional), `Toe` = toe slope (depositional) |
| `month` | – | Sampling month: `Dec`, `Jan`, `Feb` |
| `sampling_date` | date | 2023-12-15, 2024-01-18 or 2024-02-22 |
| `fraction` | – | Extractant, in the order of extraction: `H2O` (deionized water), `NaHCO3` (0.5 mol/L), `NaOH` (0.1 mol/L), `HNO3` (1 mol/L) |
| `P_mg_kg` | mg P kg⁻¹ soil | Phosphorus concentration of the fraction, mean of three analytical replicates (Table 1) |
| `P_sd_mg_kg` | mg P kg⁻¹ soil | Standard deviation of the three analytical replicates (Table 1) |
| `d18Op` | ‰ VSMOW | Oxygen isotope composition of phosphate, δ<sup>18</sup>O<sub>P</sub> |
| `d18Ow_soil_water` | ‰ VSMOW | δ<sup>18</sup>O of the soil water of the composite sample |
| `d2Hw_soil_water` | ‰ VSMOW | δ<sup>2</sup>H of the soil water of the composite sample |
| `soil_temp_C` | °C | Soil temperature at 10 cm measured in the plot at the time of sampling |

Notes

- The composite samples are in the order of the laboratory data sheet (NT before CT, top slope
  before toe slope, December to February), each with its four fractions in the order of
  extraction. `R/08_connectivity.R` resamples the composite samples by their position in this
  order, so the order of the rows must be kept for its random draws to be reproduced.
- δ<sup>18</sup>O<sub>P</sub> and the soil-water values have the precision of the data sheet
  (two or three decimals). Supplementary Table S1 prints them with two decimals.
- Soil-water δ<sup>18</sup>O, δ<sup>2</sup>H and soil temperature belong to the composite sample
  and are repeated on its four fraction rows. In one sample (CT, top slope, January) the four
  rows carry slightly different soil-water δ<sup>18</sup>O values (−13.42, −13.40, −13.38 and
  −13.38‰), as in the data sheet; each row is used for the equilibrium of its own observation,
  and their mean is used where one value per sample is needed.
- Equilibrium δ<sup>18</sup>O<sub>P</sub> and the departure from equilibrium
  (Δδ<sup>18</sup>O<sub>P</sub>) are not stored here. They are calculated by `R/01_data.R` and
  written to `output/tables/TableS1_d18Op_equilibrium.csv`.

## `precipitation_isotopes.csv`

Isotopic composition of 22 precipitation samples collected at the farm: 11 between 8 March and
17 May 2023 and 11 between 4 December 2023 and 4 June 2024.

| Column | Unit | Description |
| --- | --- | --- |
| `sampling_date` | date | Date of collection |
| `d18O` | ‰ VSMOW | δ<sup>18</sup>O of the sample |
| `d2H` | ‰ VSMOW | δ<sup>2</sup>H of the sample |

Two samples (12 and 26 April 2023) have a deuterium excess far from that of the other samples
(69.4 and −2.9‰; the other 20 lie between 4.9 and 11.0‰). `R/04_water_isotopes.R` identifies
them by that criterion and leaves them out of the meteoric water line and of Figure 2. Their
δ<sup>18</sup>O values are used for the seasonal equilibrium of Figure S1;
`output/tables/FigureS1_seasonal_equilibrium.csv` also gives the result without them.

## `pullman_2nw_daily.csv`

Daily record of the Pullman 2 NW cooperative observer station (National Weather Service
Cooperative Observer Program, station USC00456789) for 1 January 2023 to 31 December 2024
(731 days), from the monthly station summaries compiled by Washington State University. The
station reports in °F and inches; the values here are in °C (two decimals) and mm.

| Column | Unit | Description |
| --- | --- | --- |
| `date` | date | Day |
| `air_tmax_C`, `air_tmin_C` | °C | Daily maximum and minimum air temperature |
| `precip_mm` | mm | Daily precipitation |
| `soil_temp_5cm_C`, `soil_temp_10cm_C`, `soil_temp_20cm_C` | °C | Daily mean soil temperature at 5, 10 and 20 cm (2, 4 and 8 in) |

Soil temperature is missing on 32 days (1–10 January 2023 and 2–23 December 2023).
`R/02_station_climate.R` summarizes the record by sampling month, season and year
(Supplementary Table S2).
