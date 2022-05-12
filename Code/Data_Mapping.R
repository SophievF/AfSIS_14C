## Mapping data ##
## Sophie von Fromm ##
## 2022-May-05 ##

#This script contains all code to reproduce the maps that are showed in the 
#manuscript von Fromm et al. (2022)

library(tidyverse)
library(tmap)
library(sf)
library(raster)


##Load data
AfSIS_14C <- read_csv("./Data/AfSIS_data_all.csv")

##Load climate zones raster file
#Present situation
KG_p_directory <- "./Data/ClimateZones/Beck_KG_V1_present_0p0083.tif"
KG_p_raster <- raster::raster(KG_p_directory)
KG_p <- raster::extract(KG_p_raster, cbind(AfSIS_14C$Longitude,
                                           AfSIS_14C$Latitude))

#Future situation
KG_f_directory <- "./Data/ClimateZones/Beck_KG_V1_future_0p0083.tif"
KG_f_raster <- raster::raster(KG_f_directory)
KG_f <- raster::extract(KG_f_raster, cbind(AfSIS_14C$Longitude,
                                           AfSIS_14C$Latitude))

#Boundaries for Africa
data("World")
Africa_map <- World %>%
  dplyr::filter(continent == "Africa") %>% 
  st_transform(4326)

#Sub-Saharan Africa (based on definition of the United Nations geoscheme for Africa)
Africa_sp <- Africa_map %>% 
  filter(name != "Algeria",
         name != "Egypt",
         name != "Libya",
         name != "W. Sahara",
         name != "Tunisia",
         name != "Morocco",
         name != "Somalia",
         name != "Djibouti",
         name != "Sudan") %>% 
  as("Spatial")

##Crop climate zones to boundaries for SSA
#Present
KG_p_africa_crop <- raster::crop(KG_p_raster, Africa_sp)
KG_p_africa <- raster::mask(KG_p_africa_crop, Africa_sp)

#Future
KG_f_africa_crop <- raster::crop(KG_f_raster, Africa_sp)
KG_f_africa <- raster::mask(KG_f_africa_crop, Africa_sp)

#Convert AfSIS sampling locations in spatial object
AfSIS_sf <- AfSIS_14C %>% 
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

AfSIS_sf$KG_p_group <- factor(AfSIS_sf$KG_p_group,
                              labels = c("Arid", "Temperate (seasonal)",
                                         "Tropical (seasonal)", 
                                         "Temperate (humid)", 
                                         "Tropical (humid)"))

##Remove climate zones (set to NA) that are not present in SSA
#Present
KG_p_group_africa <- KG_p_africa
KG_p_group_africa[KG_p_group_africa == 0] <- NA
KG_p_group_africa[KG_p_group_africa > 17] <- NA

#Plot present situation
KG_p_group_map <- tm_shape(KG_p_group_africa) + 
  tm_raster(title = "", style = "fixed", breaks = c(1, 3, 4, 8, 14, 17),
            labels = c("Arid", "Temperate (seasonal)", "Tropical (seasonal)", 
                       "Temperate (humid)", "Tropical (humid)"),
            palette = c("#5aae61", "#9970ab", "#ffffcc", "#e7d4e8", "#d9f0d3"),
            legend.show = FALSE) +
  tm_legend(position = c(0.03, 0.05), legend.text.size = 1.4, 
            legend.width = 0.8, legend.title.size = 2, legend.height = 0.8) +
  tm_layout(main.title = "a) Current climate zones (1980-2016)", title.size = 2) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_shape(AfSIS_sf) +
  tm_dots(size = 0.5, shape = 1, jitter = 0.1) +
  tm_add_legend(type = "fill", 
                col = c("#ffffcc", "#e7d4e8", "#9970ab", "#d9f0d3", "#5aae61"),
                size = 1.5, shape = 0,
                labels = c("Arid", "Temperate (seasonal)", "Tropical (seasonal)", 
                           "Temperate (humid)",  "Tropical (humid)"),
                title = "Climate zones")
tmap_save(tm = KG_p_group_map, 
          filename = "./Figures/AfSIS_14C_FigureA1a.jpeg")

#Future situation
KG_f_group_africa <- KG_f_africa
KG_f_group_africa[KG_f_group_africa == 0] <- NA
KG_f_group_africa[KG_f_group_africa > 17] <- NA

KG_f_group_map <- tm_shape(KG_f_group_africa) + 
  tm_raster(title = "", style = "fixed", breaks = c(1, 3, 4, 8, 14, 17),
            labels = c("Tropical (humid)", "Tropical (seasonal)", "Arid",
                       "Temperate (seasonal)", "Temperate (humid)"),
            palette = c("#5aae61", "#9970ab", "#ffffcc", "#e7d4e8", "#d9f0d3"),
            legend.show = FALSE) +
  tm_layout(main.title = "b) Future climate zones (2071-2100)", title.size = 2) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_shape(AfSIS_sf) +
  tm_dots(size = 0.5, shape = 1, jitter = 0.1)
tmap_save(tm = KG_f_group_map, 
          filename = "./Figures/AfSIS_14C_FigureA1b.jpeg")
