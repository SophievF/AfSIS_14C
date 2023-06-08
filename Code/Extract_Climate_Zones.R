## Extract Köppen Geiger Climate zones for AfSIS data ##
## Sophie von Fromm ##
## 2022-04-27 ##

library(tidyverse)
library(raster)
library(tmap)
library(sf)

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
KG_f_dir <- "./Data/ClimateZones/Beck_KG_V1_future_0p0083.tif"
KG_f_raster <- raster::raster(KG_f_dir)
KG_f <- raster::extract(KG_f_raster, cbind(AfSIS_LongLat$Longitude,
                                           AfSIS_LongLat$Latitude))

# Load KG legend
KG_legend <- read_csv("./Data/ClimateZones/KG_present_legend.csv")

KG_p_legend <- KG_legend %>% 
  rename(KG_p = pro_KG_present,
         KG_p_code = pro_KG_present_short,
         KG_p_name = pro_KG_present_long)

KG_f_legend <- KG_legend %>% 
  rename(KG_f = pro_KG_present,
         KG_f_code = pro_KG_present_short,
         KG_f_name = pro_KG_present_long)

KG_p <- data.frame(KG_p) %>% 
  left_join(KG_p_legend)

KG_f <- data.frame(KG_f) %>% 
  left_join(KG_f_legend)

#Merge all three data sets
AfSIS_14C_KG <- cbind(AfSIS_LongLat, KG_p, KG_f) %>% 
  tibble() %>% 
  mutate(KG_p_group = case_when(
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

#Check which climate zones are present in africa
AfSIS_14C_KG %>% 
  group_by(KG_p_group, KG_p_code, KG_p_name) %>% 
  count()

AfSIS_14C_KG %>% 
  group_by(KG_f_group, KG_f_code, KG_f_name) %>% 
  count()

## Which Climate Zones are present across SSA
data("World")
africa_map <- World %>%
  st_transform(4326) %>% 
  dplyr::filter(continent == "Africa")

# Sub-Saharan Africa (based on definition of the United Nations geoscheme for Africa)
ssa <- africa_map %>% 
  filter(name != "Algeria",
         name != "Egypt",
         name != "Libya",
         name != "W. Sahara",
         name != "Tunisia",
         name != "Morocco",
         name != "Djibouti",
         name != "Sudan") 

ssa_sp <- ssa %>% 
  as("Spatial")

# Crop climate data to boundaries for SSA
KG_p_crop <- raster::crop(KG_p_raster, ssa_sp)
KG_p_africa <- raster::mask(KG_p_crop, ssa_sp)

KG_p_africa_sum <- KG_p_africa %>% 
  as.data.frame() %>% 
  rename(KG_p = Beck_KG_V1_present_0p0083) %>% 
  left_join(KG_p_legend) %>% 
  drop_na() %>% 
  count(KG_p, KG_p_code, KG_p_name) %>% 
  mutate(n_rel = n/sum(n)*100)
  
Climate_SSA_14C <- KG_p_africa_sum %>% 
  arrange(desc(n_rel)) %>% 
  tibble() %>% 
  dplyr::select(-KG_p, -n) %>% 
  left_join(AfSIS_14C_KG %>% 
              count(KG_p_code, KG_p_name) %>% 
              mutate(AfSIS_14C = TRUE)) %>% 
  mutate(n_rel_afsis = n/sum(n, na.rm = TRUE)*100)

write.csv(Climate_SSA_14C, row.names = FALSE,
          "./Data/Climate_SSA_14C_sum.csv")

### Extract MAT and MAP, PET
# MAT
MAT_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_01.tif"
MAT_raster <- raster(MAT_dir)
MAT <- raster::extract(MAT_raster, cbind(AfSIS_LongLat$Longitude,
                                         AfSIS_LongLat$Latitude))
# MAP
MAP_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_12.tif"
MAP_raster <- raster(MAP_dir)
MAP <- raster::extract(MAP_raster, cbind(AfSIS_LongLat$Longitude,
                                         AfSIS_LongLat$Latitude))

# PET
PET_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/global-et0_annual.tif/et0_yr/et0_yr.tif"
PET_raster <- raster(PET_dir) 
PET <- raster::extract(PET_raster, cbind(AfSIS_LongLat$Longitude,
                                         AfSIS_LongLat$Latitude))

#Merge all three data sets
AfSIS_14C_climate <- cbind(AfSIS_14C_KG, MAT, MAP, PET) %>% 
  tibble()

write.csv(AfSIS_14C_climate, row.names = FALSE,
          "./Data/AfSIS_LongLat_ClimateZones.csv")

## WRB soil types
#Jones A, Breuning-Madsen H, et al. (eds.), 2013, Soil Atlas of Africa. 
#European Commission, Publications Office of the European Union, Luxembourg. 176 pp.
sf_use_s2(FALSE)
wrb_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/Africa_soil_WRB/afticasoilmap.shp"

# Extraction soil types (WRB) from shp file####
wrb <- sf::st_read(wrb_dir)

wrb_full <- wrb %>% 
  rename(WRB_Code = SU_WRB1_PH) %>% 
  dplyr::mutate(WRB_Type = case_when(
    startsWith(WRB_Code, "AC") ~ "Acrisols",
    startsWith(WRB_Code, "AL") ~ "Alisols",
    startsWith(WRB_Code, "AN") ~ "Andosols",
    startsWith(WRB_Code, "AR") ~ "Arenosols",
    startsWith(WRB_Code, "CL") ~ "Calcisols",
    startsWith(WRB_Code, "CM") ~ "Cambisols",
    startsWith(WRB_Code, "CH") ~ "Chernozems",
    startsWith(WRB_Code, "CR") ~ "Cryosols",
    startsWith(WRB_Code, "DU") ~ "Durisols",
    startsWith(WRB_Code, "FL") ~ "Fluvisols",
    startsWith(WRB_Code, "FR") ~ "Ferralsols",
    startsWith(WRB_Code, "GL") ~ "Gleysols", 
    startsWith(WRB_Code, "HS") ~ "Histosols",
    startsWith(WRB_Code, "KS") ~ "Kastanozems",
    startsWith(WRB_Code, "LP") ~ "Leptosols",
    startsWith(WRB_Code, "LX") ~ "Lixisols",
    startsWith(WRB_Code, "LV") ~ "Luvisols",
    startsWith(WRB_Code, "NT") ~ "Nitisols",
    startsWith(WRB_Code, "PH") ~ "Phaeozems",
    startsWith(WRB_Code, "PL") ~ "Planosols",
    startsWith(WRB_Code, "PT") ~ "Plinthosols",
    startsWith(WRB_Code, "PZ") ~ "Podzols",
    startsWith(WRB_Code, "RG") ~ "Regosols",
    startsWith(WRB_Code, "SC") ~ "Solonchaks",
    startsWith(WRB_Code, "SN") ~ "Solonetz",
    startsWith(WRB_Code, "ST") ~ "Stagnosols",
    startsWith(WRB_Code, "TC") ~ "Technosols",
    startsWith(WRB_Code, "UM") ~ "Umbrisols", 
    startsWith(WRB_Code, "VR") ~ "Vertisols",
    startsWith(WRB_Code, "GY") ~ "Gypsisols",
    startsWith(WRB_Code, "WR") ~ "Water body"
  ))

# merging afsis data with wrb
afsis_sf <- sf::st_as_sf(AfSIS_LongLat, coords = c("Longitude", "Latitude"), crs = 4326)

afsis_wrb_sf <- sf::st_join(afsis_sf, wrb_full, left = TRUE)

wrb_afsis <- st_set_geometry(afsis_wrb_sf, NULL)

wrb_afsis %>% 
  count(WRB_Type)

#reduce area to ssa

wrb_ssa <- wrb_full[ssa,]

wrb_ssa_sum <- wrb_ssa %>% 
  data.frame() %>% 
  drop_na(WRB_Type) %>% 
  count(WRB_Type)

wrb_ssa_14c <- wrb_ssa_sum %>% 
  filter(WRB_Type != "Water body") %>% 
  mutate(n_rel = n/sum(n)*100) %>% 
  dplyr::select(-n) %>% 
  left_join(wrb_afsis %>% 
              count(WRB_Type ) %>% 
              mutate(AfSIS_14C = TRUE)) %>% 
  mutate(n_rel_afsis = n/sum(n, na.rm = TRUE)*100) %>% 
  arrange(desc(n_rel))
write.csv(wrb_ssa_14c, row.names = FALSE,
          "./Data/WRB_SSA_14C_sum.csv")

