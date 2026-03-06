# BioMeteoLab NEE Modelling

## Data Merge script (`Data_Merge.R`)

Reads multi-sheet collar data (16 sensors) and 10-min logger data, aligns timestamps to a 10-min grid, reshapes collars to wide format, and merges with logger.  
**Output:** `data/Merged_data.csv`

### Source data

Both datasets are available on the **LTOM drive**, folder **NEE_Modelling_Data**, subfolder **Initial Data**. The files there are slightly edited versions of the originals from Sharepoint. Use these.

### How to run

1. Create a folder named `data` in the project root (same level as `Data_Merge.R`).
2. Copy the two Excel files (Iisaku collars and Iisaku logger) from that location into the `data` folder.
3. Run `Data_Merge.R`.

## NEE Random Forest by collar (`NEE_RF_Collar_Analysis.R`)

Fits a separate Random Forest model for each collar (1–16): response is that collar’s NEE (`NEE.CO2_collars_i`), predictors are that collar’s local variables (PAR, temp, moisture) plus shared logger variables (soil/air temperature, VWC, PA_uS, rain, hour, doy). Uses **ranger** (500 trees, permutation importance), 50/50 train–test split. No CSV or RDS is written; metrics and variable importance are printed to the console. Plots are saved only in `plots/` (R² by collar, RMSE by collar, heatmap of variable importance with top‑3 factors per collar marked).

**Requires:** `tidyverse`, `lubridate`, `ranger`, `ggplot2`.

**How to run:** Ensure `data/Merged_data.csv` exists (run `Data_Merge.R` first). Run `NEE_RF_Collar_Analysis.R`; it loads `merged_data` from that CSV if not already in the R environment. The script creates `plots/` if missing and saves `plots/rf_r2_by_collar.png`, `plots/rf_importance_heatmap_by_collar.png`, `plots/rf_rmse_by_collar.png`.
