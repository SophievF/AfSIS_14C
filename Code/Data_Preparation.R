## Sophie von Fromm
## 2022-May-04

library(tidyverse)

## Data preparation - bring all data together
# Radiocarbon data (Load_14C_Data.R, MeanC_age.R)
# XRPD data (XRPD_14C_MineralFitting.R)
# Soil and field data (csv file)
# Climate zones (Extract_Climate_Zones.R)
# GPP (Gross_Primary_Productivity.R)

## This is a subset of the AfSIS Phase I dataset
# More information can be found in von Fromm et al. (2021), SOIL (https://doi.org/10.5194/soil-7-305-2021)
AfSIS_soil_data <- read_csv("./Data/AfSIS_Soil_Data.csv")

## Load all the other dataset
# Radiocarbon
AfSIS_14C <- read_csv("./Data/df_14C_TT.csv")

# XRPD data
AfSIS_XRPD <- read_csv("./Data/AfSIS_Mineral_fits.csv")

# Climate zones
AfSIS_climate <- read_csv("./Data/AfSIS_LongLat_ClimateZones.csv")

# GPP
AfSIS_GPP <- read_csv("./Data/AfSIS_GPP.csv") %>% 
  dplyr::rename(SSN = SiteID) %>% 
  #Coordinates have one digit less than orginal data
  dplyr::select(SSN, GPP)

AfSIS_all <- AfSIS_soil_data %>% 
  left_join(AfSIS_GPP) %>% 
  left_join(AfSIS_climate) %>% 
  left_join(AfSIS_14C) %>% 
  left_join(AfSIS_XRPD)

## Prepare data for analysis

AfSIS_data <- AfSIS_all %>% 
  # Erosion: was measured at 4 sub-plots (Erosion1-Erosion4); need to calculate for 
  # plot level; 0 - no erosion, 1 - sheet, 2 - rill 3 - gully
  # Different cut-offs (1,2,3) result in almost the same classification
  rowwise() %>% 
  mutate(Erosion = case_when(
    Erosion1 + Erosion2 + Erosion3 + Erosion4 > 2 ~ "yes",
    TRUE ~ "no" 
  )) %>% 
  ungroup() %>% 
  # remove sub-plot erosion data
  dplyr::select(-c(Erosion1:Erosion4)) %>% 
  # convert cultivation in factor
  mutate(PlotCultMgd = case_when(
    PlotCultMgd == 0 ~ "no",
    TRUE ~ "yes"
  )) %>% 
  # create mineral groups
  rowwise() %>% 
  mutate(Clay_1_1 = Kaolinite + Dickite + Halloysite,
         Clay_2_1 = `Smectite (ML)` + `Smectite (Di)` + Illite + Vermiculite ,
         Amorphous = Glass + Ferrihydrite,
         Fehydroxide = Goethite + Hematite,
         Feldspars = `K-feldspar` + Plagioclase) %>%
  ungroup() %>% 
  # remove other minerals
  dplyr::select(-c(Kaolinite:Background)) %>% 
  dplyr::mutate(Mox = Alox + (1/2 * Feox))

# Check 14C data for anomalies
plotly::ggplotly(
  AfSIS_data %>% 
    ggplot(aes(y = Delta14C, x = CORG, color = LandCover, group = SSN)) +
    geom_point() +
    facet_wrap(~Depth)
)

# samples that were measured with gas inlet have no 13C data
plotly::ggplotly(
  AfSIS_data %>% 
    ggplot(aes(y = Delta14C, x = Delta13C, color = LandCover, group = SSN)) +
    geom_point() +
    facet_wrap(~Depth)
)

# Two samples (icr042266 and icr042245)  show unusally high 13C values 
# These samples probably contain inorganic C which is supported by the low 14C values
# Because of this we exclude the plots from which the samples are

AfSIS_data <- AfSIS_data %>% 
  filter(SSN != "icr042245",
         SSN != "icr042246",
         SSN != "icr042265",
         SSN != "icr042266")

write.csv(AfSIS_data, row.names = FALSE,
          "./Data/AfSIS_data_all.csv")


