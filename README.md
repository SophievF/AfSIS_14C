# AfSIS_14C
Radiocarbon analaysis of a subset of the AfSIS dataset (related publications: von Fromm et al. (2021): https://soil.copernicus.org/articles/7/305/2021/ and Vagen et al. (2021): https://data.worldagroforestry.org/dataset.xhtml?persistentId=doi:10.34725/DVN/66BFOB

Radiocarbon data will also be made publicily availabe via the Internation Soil Radiocarbon Database: https://github.com/International-Soil-Radiocarbon-Database/ISRaD

Authors: Sophie von Fromm and Benjamin Butler

Date: November 2023 

This repository contains all the code to reproduce the analysis and all figures in the publication von Fromm et al. (2023).

Only the file AfSIS_data_all, and the R scripts Data_Mapping.R, Data_Analysis.R, DataRepresentative_Analysis.R, SpectralData_Analysis.R are needed to reproduce the analysis and figures in the manuscript. All other files either contain raw data and/or are part of the data preparation and can be generated with the corresponding R scripts.

The folder Code contains all the R code, the folder Data contains all the data needed to run the R Scripts and the folder Figure contains all the figures and tables that are produced with the R code.

Folder Data:
  - AfSIS_data_all.csv: Merged data file with all data to reproduce analysis and figures. Output file from Data_Preparation.R
  - AfSIS_LongLat.csv: Longitute and latitute data for the AfSIS samples
  - AfSIS_XRPD_LongLat.csv: Longitute and latitute data for the X-ray powder diffraction (XRPD) samples
  - AfSIS_Soil_data.csv: Chemical and field data for the AfSIS samples. https://doi.org/10.34725/DVN/66BFOB
  - Folder ClimateZones: Geospatial data for present and future climate zones from Beck et al. 2018, Sci Data, https://doi.org/10.1038/sdata.2018.214
  - Folder GPP: Gross Primary Productivity data from the Fluxcom network (https://www.fluxcom.org/CF-Download/) for 01/2001 to 12/2019 for the AfSIS sampling locations 
  - Folder Radiocarbon: Raw radiocarbon data from the 14C analysis lab at MPI-BGC, Jena
  - Folder XRPD_XY_Files: X-ray powder diffraction patterns stored as XY files for each sample
  - dd_14C_TT.csv: Outpfile from MeanC_age.R
  - AfSIS_Mineral_fits.csv: Output file from XRPD_14C_MineralFitting.R
  - AfSIS_GPP.csv: Output file from Gross_Primary_Productivity.R
  - AfSIS_LongLat_ClimateZones.csv: Output file from Extract_Climate_Zones.R

Folder Code:
  - Data_Preparation.R: Script to merge all data files. Outpfile: AfSIS_data_all.csv
  - DataRepresentative_Analysis.R. Script to perform the representative analysis of the different AfSIS datasets (see suplement in manuscript)
  - Data_Analysis.R: Script to reproduce analysis and figures
  - Data_Analysis_LinearRegression.R: Script to reproduce linear mixed-effect models analysis (including related figures and tables)
  - Data_AnaLysis_RandomForest.R: Script to reproduce random forest model anaylsis including partial dependence plots
  - Data_Mapping.R: Script to reproduce all maps
  - Load14C_Data.R: Function to load raw 14C data
  - MeanC_age.R: Script to calculate mean C age / turnover time. Output file: dd_14C_TT.csv
  - XRPD_14C_MineralFitting.R: Script to quantify mineral content in each sample. Output file: AfSIS_Mineral_fits.csv
  - Gross_Primary_Productivity.R: Script to calculate gross primary productivity. Output file: AfSIS_GPP.csv
  - Extract_Climate_Zones: Script to extract climate zones from global product for each sampling location. Output file: AfSIS_LongLat_ClimateZones.csv
