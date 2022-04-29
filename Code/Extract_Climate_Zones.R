## Extract Köppen Geiger Climate zones for AfSIS data ##
## Sophie von Fromm ##
## 2022-04-27 ##

library(tidyverse)
library(raster)

# KG data is from Beck et al. 2018
# https://doi.org/10.1038/sdata.2018.214
# The data can be accessed via: https://doi.org/10.6084/m9.figshare.6396959

# p - present: 1980-2016
# f - future: 2071-2100 (RCP8.5)

# Load AfSiS data
AfSIS_LongLat <- read_csv("./Data/AfSIS_LongLat.csv")

# Load present climate zones
KG_p_dir <- "./Data/ClimateZones/Beck_KG_V1_present_0p0083.tif"
KG_p_raster <- raster::raster(KG_p_dir)
KG_p <- raster::extract(KG_p_raster, cbind(AfSIS_LongLat$Longitude,
                                           AfSIS_LongLat$Latitude))

# Load projected climate zones
KG_f_dir <- ".Data/ClimateZones/Beck_KG_V1_future_0p0083.tif"
KG_f_raster <- raster::raster(KG_f_dir)
KG_f <- raster::extract(KG_f_raster, cbind(AfSIS_LongLat$Longitude,
                                           AfSIS_LongLat$Latitude))

#Merge all three data sets
AfSIS_14C_KG <- cbind(AfSIS_LongLat, KG_p, KG_f) %>% 
  tibble()

#Check which climate zones are present in africa
AfSIS_14C_KG %>% 
  group_by(KG_p) %>% 
  count()

AfSIS_14C_KG %>% 
  group_by(KG_f) %>% 
  count()

#Convert numbers in actual climate zones based on legend.txt
#Create new group (KG_p/f_group) that summarizes main climate zones
AfSIS_14C_KG <- AfSIS_14C_KG %>% 
  mutate(KG_p_code = case_when(
    KG_p == 1 ~ "Af",
    KG_p == 2 ~ "Am",
    KG_p == 3 ~ "Aw",
    KG_p == 4 ~ "BWh",
    KG_p == 6 ~ "BSh",
    KG_p == 7 ~ "BSk",
    KG_p == 9 ~ "Csb",
    KG_p == 11 ~ "Cwa",
    KG_p == 12 ~ "Cwb",
    KG_p == 14 ~ "Cfa",
    KG_p == 15 ~ "Cfb"
  ),
  KG_f_code = case_when(
    KG_f == 1 ~ "Af",
    KG_f == 2 ~ "Am",
    KG_f == 3 ~ "Aw",
    KG_f == 4 ~ "BWh",
    KG_f == 6 ~ "BSh",
    KG_f == 11 ~ "Cwa",
    KG_f == 12 ~ "Cwb",
    KG_f == 15 ~ "Cfb"
  ),
  KG_p_name = case_when(
    KG_p == 1 ~ "Tropical, rainforest",
    KG_p == 2 ~ "Tropical, monsoon",
    KG_p == 3 ~ "Tropical, savannah",
    KG_p == 4 ~ "Arid, desert, hot",
    KG_p == 6 ~ "Arid, steppe, hot",
    KG_p == 7 ~ "Arid, steppe, cold",
    KG_p == 9 ~ "Temperate, dry summer, warm summer",
    KG_p == 11 ~ "Temperate, dry winter, hot summer",
    KG_p == 12 ~ "Temperate, dry winter, warm summer",
    KG_p == 14 ~ "Temperate, no dry season, hot summer",
    KG_p == 15 ~ "Temperate, no dry season, warm summer"
  ),
  KG_f_name = case_when(
    KG_f == 1 ~ "Tropical, rainforest",
    KG_f == 2 ~ "Tropical, monsoon",
    KG_f == 3 ~ "Tropical, savannah",
    KG_f == 4 ~ "Arid, desert, hot",
    KG_f == 6 ~ "Arid, steppe, hot",
    KG_f == 11 ~ "Temperate, dry winter, hot summer",
    KG_f == 12 ~ "Temperate, dry winter, warm summer",
    KG_f == 15 ~ "Temperate, no dry season, warm summer"
  ),
  KG_p_group = case_when(
    KG_p == 1 ~ "Tropical (humid)",
    KG_p == 2 ~ "Tropical (humid)",
    KG_p == 3 ~ "Tropical (seasonal)",
    KG_p == 4 ~ "Arid",
    KG_p == 6 ~ "Arid",
    KG_p == 7 ~ "Arid",
    KG_p == 9 ~ "Temperate (seasonal)",
    KG_p == 11 ~ "Temperate (seasonal)",
    KG_p == 12 ~ "Temperate (seasonal)",
    KG_p == 14 ~ "Temperate (humid)",
    KG_p == 15 ~ "Temperate (humid)"
  ),
  KG_f_group = case_when(
    KG_f == 1 ~ "Tropical (humid)",
    KG_f == 2 ~ "Tropical (humid)",
    KG_f == 3 ~ "Tropical (seasonal)",
    KG_f == 4 ~ "Arid",
    KG_f == 6 ~ "Arid",
    KG_f == 11 ~ "Temperate (seasonal)",
    KG_f == 12 ~ "Temperate (seasonal)",
    KG_f == 15 ~ "Temperate (humid)"
  )) %>% 
  dplyr::select(-KG_p, -KG_f)

write.csv(AfSIS_14C_XRPD_Global_KG, row.names = FALSE,
          "./Data/AfSIS_LongLat_ClimateZones.csv")
