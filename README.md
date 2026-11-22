# Data

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
| `d18Op` | ‰ VSMOW | Oxygen isotope composition of phosphate, δ<sup>18</sup>O<sub>P</sub> (Supplementary Table S1) |
| `d18Ow_soil_water` | ‰ VSMOW | δ<sup>18</sup>O of the soil water of the composite sample (Supplementary Table S1) |
| `soil_temp_C` | °C | Soil temperature at 10 cm measured in the plot at the time of sampling |

Soil-water δ<sup>18</sup>O and soil temperature belong to the composite sample and are repeated
on its four fraction rows, as in Supplementary Table S1.

Equilibrium δ<sup>18</sup>O<sub>P</sub> and the departure from equilibrium
(Δδ<sup>18</sup>O<sub>P</sub>) are not stored here. They are calculated by `R/01_data.R` and
written to `output/tables/TableS1_d18Op_equilibrium.csv`.

## `soil_temperature_summary.csv`

Daily mean soil temperature at the Pullman 2 NW cooperative observer station (USC00456789),
1 January 2023 to 31 December 2024, summarized by period and depth (Supplementary Table S2).

| Column | Unit | Description |
| --- | --- | --- |
| `period` | – | Month, season or year over which the daily means are summarized |
| `depth_cm` | cm | Sensor depth: 5, 10 or 20 cm (2, 4 and 8 in) |
| `min_C`, `mean_C`, `max_C` | °C | Minimum, mean and maximum of the daily mean soil temperature in the period |
| `days_with_data` | days | Days in the period with a soil-temperature record |
| `days_total` | days | Days in the period |

The three sampling-month rows give the temperature ranges of the equilibrium bands in Figure 3;
the row for the sampling months combined is used in Supplementary Text S1. The daily station
record itself is not included; it is available from the NOAA National Centers for Environmental
Information (Global Historical Climatology Network – Daily).
