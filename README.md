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
