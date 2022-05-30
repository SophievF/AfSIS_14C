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

####Future situation
###Differences between present and future KG zones
##Merge KG sub-groups in 5 main groups 
#1-arid, 2-temperate (dry), 3-temperate (humid), 4-tropical (dry), 5-tropical (humid)
#values >0 & <= 2 ~ 5; >2 & <= 3 ~ 4; >3 & <= 7 ~ 1; >7 & <= 13 ~ 2; < 13 & <= 16 ~ 3

recal <- c(0,2,5, 2,3,4, 3,7,1, 7,13,2, 13,16,3)
recal_mat <- matrix(recal, ncol = 3, byrow = TRUE)

KG_p_africa_grouped <- reclassify(KG_p_group_africa, recal_mat)

##Areas that become arid ~ 1; (only keep arid values)
recal_NA_arid <- c(0,2,NA, 2,3,NA, 3,7,1, 7,13,NA, 13,16,NA)
recal_mat_NA_arid <- matrix(recal_NA_arid, ncol = 3, byrow = TRUE)

#Only keep areas that are arid in the future (NA values are removed)
KG_f_africa_new_arid <- reclassify(KG_f_group_africa, recal_mat_NA_arid)
plot(KG_f_africa_new_arid)

#Areas that become arid are highlighted (present - future)
KG_p_f_arid <- overlay(KG_p_africa_grouped,
                       KG_f_africa_new_arid,
                       fun = function(r1, r2){return(r1-r2)})
plot(KG_p_f_arid)

##New groups/values
#arid ~ 1 - arid ~ 1 --> 0 (no change)
#temperate (dry) ~ 2 - arid ~ 1 --> 1
#temperate (humid) ~ 3 - arid ~ 1 --> 2
#tropical (dry) ~ 4 - arid ~ 1 --> 3
#tropical (humid) ~ 5 - arid ~ 1 --> 4 (doesn't exist)

map_arid <- tm_shape(KG_p_f_arid) + 
  tm_raster(title = "Climate becomes arid", style = "cat",
            labels = c("no change",
                       "temperate (dry season) -> arid",
                       "temperate (humid) -> arid",
                       "tropical (dry season) -> arid"),
            palette = c("darkgrey",
                        "#fff7bc",
                        "#fec44f",
                        "#d95f0e")) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_legend(legend.position = c("left", "bottom"),
            legend.width = 1)

##Areas that become temperate (dry) ~ 2; (only keep temperate (dry))
recal_NA_tempdry <- c(0,2,NA, 2,3,NA, 3,7,NA, 7,13,2, 13,16,NA)
recal_mat_NA_tempdry <- matrix(recal_NA_tempdry, ncol = 3, byrow = TRUE)

KG_f_africa_new_tempdry <- reclassify(KG_f_group_africa, recal_mat_NA_tempdry)
plot(KG_f_africa_new_tempdry)

KG_p_f_tempdry <- overlay(KG_p_africa_grouped,
                          KG_f_africa_new_tempdry,
                          fun = function(r1,r2){return(r1-r2)})
plot(KG_p_f_tempdry)

##New groups/values
#temperate (dry) ~ 2 ~ 1 - temperate (dry) ~ 2 --> 0 (no change) 
#arid ~ 1 - temperate (dry) ~ 2 --> -1
#temperate (humid) ~ 3 - temperate (dry) ~ 2 --> 1
#tropical (dry) ~ 4 - temperate (dry) ~ 2 --> 2 (doesn't exist)
#tropical (humid) ~ 5 - temperate (dry) ~ 2 --> 3 (doesn't exist)

map_tempdry <- tm_shape(KG_p_f_tempdry) + 
  tm_raster(title = "Climate becomes temperate (dry season)", style = "cat",
            labels = c("arid -> temperate (dry season)",
                       "no change",
                       "temperate (humid) -> temperate (dry season)"),
            palette = c("#fde0dd",
                        "darkgrey",
                        "#fa9fb5")) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_legend(legend.position = c("left", "bottom"),
            legend.width = 1)

##Areas that become temperate (humid) ~ 3; (only keep temperate (humid))
recal_NA_tempwet <- c(0,2,NA, 2,3,NA, 3,7,NA, 7,13,NA, 13,16,3)
recal_mat_NA_tempwet <- matrix(recal_NA_tempwet, ncol = 3, byrow = TRUE)

KG_f_africa_new_tempwet <- reclassify(KG_f_group_africa, recal_mat_NA_tempwet)
plot(KG_f_africa_new_tempwet)

KG_p_f_tempwet <- overlay(KG_p_africa_grouped,
                          KG_f_africa_new_tempwet,
                          fun = function(r1,r2){return(r1-r2)})
plot(KG_p_f_tempwet)

##New groups/values
#temperate (humid) ~ 3 - temperate (humid) ~ 3 --> 0 (no change) 
#arid ~ 1 - temperate (humid) ~ 3 --> -2 (doesn't exist)
#temperate (dry) ~ 2 - temperate (humid) ~ 3 --> -1
#tropical (dry) ~ 4 - temperate (humid) ~ 3 --> 1 (doesn't exist)
#tropical (humid) ~ 5 - temperate (humid) ~ 3 --> 2 (doesn't exist)

map_tempwet <- tm_shape(KG_p_f_tempwet) + 
  tm_raster(title = "Climate becomes temperate (humid)", style = "cat",
            labels = c("temperate (dry season) -> temperate (humid)",
                       "no change"),
            palette = c("#a8ddb5",
                        "darkgrey")) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_legend(legend.position = c("left", "bottom"),
            legend.width = 1)

##Areas that become tropical (dry) ~ 4; (only keep tropical (dry))
recal_NA_tropdry <- c(0,2,NA, 2,3,4, 3,7,NA, 7,13,NA, 13,16,NA)
recal_mat_NA_tropdry <- matrix(recal_NA_tropdry, ncol = 3, byrow = TRUE)

KG_f_africa_new_tropdry <- reclassify(KG_f_group_africa, recal_mat_NA_tropdry)
plot(KG_f_africa_new_tropdry)

KG_p_f_tropdry <- overlay(KG_p_africa_grouped,
                          KG_f_africa_new_tropdry,
                          fun = function(r1,r2){return(r1-r2)})
plot(KG_p_f_tropdry)

##New groups/values
#tropical (dry) ~ 4 - tropical (dry) ~ 4 --> 0 (no change) 
#arid ~ 1 - tropical (dry) ~ 4 --> -3 
#temperate (dry) ~ 2 - tropical (dry) ~ 4 --> -2
#temperate (humid) ~ 3 - tropical (dry) ~ 4 --> -1 
#tropical (humid) ~ 5 - tropical (dry) ~ 4 --> 1 

map_tropdry <- tm_shape(KG_p_f_tropdry) + 
  tm_raster(title = "Becomes tropical (dry season)", style = "cat",
            labels = c("arid -> tropical (dry season)",
                       "temperate (dry season) -> tropical (dry season)",
                       "temperate (humid) -> tropical (dry season)",
                       "no change",
                       "tropical (humid) -> tropical (dry season)"),
            palette = c("#edf8e9",
                        "#bae4b3",
                        "#74c476",
                        "darkgrey",
                        "#238b45")) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_legend(legend.position = c("left", "bottom"),
            legend.width = 1)

##Areas that become tropical (humid) ~ 5; (only keep tropical (humid))
recal_NA_tropwet <- c(0,2,5, 2,3,NA, 3,7,NA, 7,13,NA, 13,16,NA)
recal_mat_NA_tropwet <- matrix(recal_NA_tropwet, ncol = 3, byrow = TRUE)

KG_f_africa_new_tropwet <- reclassify(KG_f_group_africa, recal_mat_NA_tropwet)
plot(KG_f_africa_new_tropwet)

KG_p_f_tropwet <- overlay(KG_p_africa_grouped,
                          KG_f_africa_new_tropwet,
                          fun = function(r1,r2){return(r1-r2)})
plot(KG_p_f_tropwet)

##New groups/values
#tropical (humid) ~ 5 - tropical (humid) ~ 5 --> 0 (no change) 
#arid ~ 1 - tropical (humid) ~ 5 --> -4 (doesn't exist)
#temperate (dry) ~ 2 - tropical (humid) ~ 5 --> -3
#temperate (humid) ~ 3 - tropical (humid) ~ 5 --> -2 
#tropical (dry) ~ 4 - tropical (humid) ~ 5 --> -1 

map_tropwet <- tm_shape(KG_p_f_tropwet) + 
  tm_raster(title = "Becomes tropical (humid)", style = "cat",
            labels = c("temperate (dry season) -> tropical (humid)",
                       "temperate (humid) -> tropical (humid)",
                       "tropical (dry season) -> tropical (humid)",
                       "no change"),
            palette = c("#deebf7",
                        "#9ecae1",
                        "#3182bd",
                        "darkgrey")) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "#646464", lwd = 0.5) +
  tm_legend(legend.position = c("left", "bottom"),
            legend.width = 1)

###Create map with 'highest potential to stabilize C': temperate dry & tropical dry
##Figure 4
recal_NA_all <- c(0,1,NA, 1,2,1, 2,3,NA, 3,4,1, 4,5,NA)
recal_mat_NA_all <- matrix(recal_NA_all, ncol = 3, byrow = TRUE)
KG_p_africa_grouped_new <- reclassify(KG_p_africa_grouped, recal_mat_NA_all)
plot(KG_p_africa_grouped_new)

##Areas with high potential to stabilize C that become arid
#Based on all areas that become arid
#arid ~ 1 - arid ~ 1 --> 0 (no change) --> NA
#temperate (dry) ~ 2 - arid ~ 1 --> 1 --> 1
#temperate (humid) ~ 3 - arid ~ 1 --> 2 --> NA
#tropical (dry) ~ 4 - arid ~ 1 --> 3 --> 1
#tropical (humid) ~ 5 - arid ~ 1 --> 4 (doesn't exist)

recal_NA_arid <- c(-0.1,0,NA, 0,1,1, 1,2,NA, 2,3,1)
recal_mat_NA_arid <- matrix(recal_NA_arid, ncol = 3, byrow = TRUE)
KG_p_f_arid_NA <- reclassify(KG_p_f_arid, recal_mat_NA_arid)
plot(KG_p_f_arid_NA)

##Areas with high potential to stabilize C that become temperate (dry)
#Based on all areas that become temperate (dry)
#Doesn't exist

##Areas with high potential to stabilize C that become temperate (humid)
#Based on all areas that become temperate (humid)
#temperate (humid) ~ 3 - temperate (humid) ~ 3 --> 0 (no change) --> NA 
#arid ~ 1 - temperate (humid) ~ 3 --> -2 (doesn't exist)
#temperate (dry) ~ 2 - temperate (humid) ~ 3 --> -1 --> 1
#tropical (dry) ~ 4 - temperate (humid) ~ 3 --> 1 (doesn't exist)
#tropical (humid) ~ 5 - temperate (humid) ~ 3 --> 2 (doesn't exist)
recal_NA_tempwet <- c(-2,-1,1, -1,2,NA)
recal_mat_NA_tempwet <- matrix(recal_NA_tempwet, ncol = 3, byrow = TRUE)
KG_p_f_tempwet_NA <- reclassify(KG_p_f_tempwet, recal_mat_NA_tempwet)
plot(KG_p_f_tempwet_NA)

##Areas with high potential to stabilize C that become tropical (dry)
#Based on all areas that become tropical (dry)
#tropical (dry) ~ 4 - tropical (dry) ~ 4 --> 0 (no change) --> NA
#arid ~ 1 - tropical (dry) ~ 4 --> -3 --> NA
#temperate (dry) ~ 2 - tropical (dry) ~ 4 --> -2 --> 1
#temperate (humid) ~ 3 - tropical (dry) ~ 4 --> -1 --> NA
#tropical (humid) ~ 5 - tropical (dry) ~ 4 --> 1 --> NA

recal_NA_tropdry <- c(-4,-3,NA, -3,-2,1, -2,4,NA)
recal_mat_NA_tropdry <- matrix(recal_NA_tropdry, ncol = 3, byrow = TRUE)
KG_p_f_tropdry_NA <- reclassify(KG_p_f_tropdry, recal_mat_NA_tropdry)
plot(KG_p_f_tropdry_NA)

##Areas with high potential to stabilize C that become tropical (humid)
#Based on all areas that become tropical (humid)
#tropical (humid) ~ 5 - tropical (humid) ~ 5 --> 0 (no change)  --> NA
#arid ~ 1 - tropical (humid) ~ 5 --> -4 (doesn't exist) 
#temperate (dry) ~ 2 - tropical (humid) ~ 5 --> -3 --> 1
#temperate (humid) ~ 3 - tropical (humid) ~ 5 --> -2 --> NA
#tropical (dry) ~ 4 - tropical (humid) ~ 5 --> -1 --> 1

recal_NA_tropwet <- c(-4,-3,1, -3,-2,NA, -2,-1,1, -1,0,NA)
recal_mat_NA_tropwet <- matrix(recal_NA_tropwet, ncol = 3, byrow = TRUE)
KG_p_f_tropwet_NA <- reclassify(KG_p_f_tropwet, recal_mat_NA_tropwet)
plot(KG_p_f_tropwet_NA)

##All other changes (inside and outside of region with highest potential to stabilize C by minerals)
#Based on all areas that become arid
#arid ~ 1 - arid ~ 1 --> 0 (no change) --> NA
#temperate (dry) ~ 2 - arid ~ 1 --> 1 --> NA
#temperate (humid) ~ 3 - arid ~ 1 --> 2 --> 1
#tropical (dry) ~ 4 - arid ~ 1 --> 3 --> NA
#tropical (humid) ~ 5 - arid ~ 1 --> 4 (doesn't exist)

recal_NA_arid_oth <- c(-0.1,0,NA, 0,1,NA, 1,2,1, 2,3,NA)
recal_mat_NA_arid_oth <- matrix(recal_NA_arid_oth, ncol = 3, byrow = TRUE)
KG_p_f_arid_NA_oth <- reclassify(KG_p_f_arid, recal_mat_NA_arid_oth)
plot(KG_p_f_arid_NA_oth)

##Areas with high potential to stabilize C that become temperate (dry)
#Based on all areas that become temperate (dry)
#arid ~ 1 - temperate (dry) ~ 2 --> -1 --> 1
#temperate (dry) ~ 2 - temperate (dry) ~ 2 --> 0 (no change) --> NA
#temperate (humid) ~ 3 - temperate (dry) ~ 2 --> 1 --> 1
#tropical (dry) ~ 4 - temperate (dry) ~ 2 --> 3 --> (doesn't exist)
#tropical (humid) ~ 5 - temperate (dry) ~ 2 --> 4 (doesn't exist)

recal_NA_tempdry_oth <- c(-2,-1,1, -0.1,0,NA, 0,1,1)
recal_mat_NA_tempdry_oth <- matrix(recal_NA_tempdry_oth, ncol = 3, byrow = TRUE)
KG_p_f_tempdry_NA_oth <- reclassify(KG_p_f_tempdry, recal_mat_NA_tempdry_oth)
plot(KG_p_f_tempdry_NA_oth)

#Based on all areas that become temperate (humid)
#temperate (humid) ~ 3 - temperate (humid) ~ 3 --> 0 (no change) --> NA 
#arid ~ 1 - temperate (humid) ~ 3 --> -2 (doesn't exist)
#temperate (dry) ~ 2 - temperate (humid) ~ 3 --> -1 --> NA
#tropical (dry) ~ 4 - temperate (humid) ~ 3 --> 1 (doesn't exist)
#tropical (humid) ~ 5 - temperate (humid) ~ 3 --> 2 (doesn't exist)
##no other changes for this climate zone

#Based on all areas that become tropical (dry)
#tropical (dry) ~ 4 - tropical (dry) ~ 4 --> 0 (no change) --> NA
#arid ~ 1 - tropical (dry) ~ 4 --> -3 --> 1
#temperate (dry) ~ 2 - tropical (dry) ~ 4 --> -2 --> NA
#temperate (humid) ~ 3 - tropical (dry) ~ 4 --> -1 --> 1
#tropical (humid) ~ 5 - tropical (dry) ~ 4 --> 1 --> 1

recal_NA_tropdry_oth <- c(-4,-3,1, -3,-2,NA, -2,-1,1, -0.1,0,NA, 0,1,1)
recal_mat_NA_tropdry_oth <- matrix(recal_NA_tropdry_oth, ncol = 3, byrow = TRUE)
KG_p_f_tropdry_NA_oth <- reclassify(KG_p_f_tropdry, recal_mat_NA_tropdry_oth)
plot(KG_p_f_tropdry_NA_oth)

#Based on all areas that become tropical (humid)
#tropical (humid) ~ 5 - tropical (humid) ~ 5 --> 0 (no change)  --> NA
#arid ~ 1 - tropical (humid) ~ 5 --> -4 (doesn't exist) 
#temperate (dry) ~ 2 - tropical (humid) ~ 5 --> -3 --> NA
#temperate (humid) ~ 3 - tropical (humid) ~ 5 --> -2 --> 1
#tropical (dry) ~ 4 - tropical (humid) ~ 5 --> -1 --> NA

recal_NA_tropwet_oth <- c(-4,-3,NA, -3,-2,1, -2,-1,NA, -1,0,NA)
recal_mat_NA_tropwet_oth <- matrix(recal_NA_tropwet_oth, ncol = 3, byrow = TRUE)
KG_p_f_tropwet_NA_oth <- reclassify(KG_p_f_tropwet, recal_mat_NA_tropwet_oth)
plot(KG_p_f_tropwet_NA_oth)

#Merge all rasters that contain 'other' changes
KG_p_f_other <- merge(KG_p_f_arid_NA_oth, KG_p_f_tempdry_NA_oth,
                      KG_p_f_tropdry_NA_oth, KG_p_f_tropwet_NA_oth)

#Merge all rasters that contain 'other' changes
KG_p_f_other <- merge(KG_p_f_arid_NA_oth, KG_p_f_tempdry_NA_oth,
                      KG_p_f_tropdry_NA_oth, KG_p_f_tropwet_NA_oth)

#Create polygon for temperate and tropical dry season regions
KG_p_africa_poly <- terra::as.polygons(as(KG_p_africa_grouped_new, "SpatRaster"))

crs(KG_p_africa_poly) <- "+proj=longlat +datum=WGS84"

plot(KG_p_africa_poly)

terra::writeVector(KG_p_africa_poly, "./Data/ClimateZones/KG_p_africa_poly.shp")

#Load polygon
sf_use_s2(FALSE)
KG_p_africa_poly_shp <- sf::read_sf("./Data/ClimateZones/KG_p_africa_poly.shp")

#Plot predicted changes
KG_p_f_NA <- tm_shape(Africa_sp, projection = 4326) +
  tm_borders(lwd = 0.5) +
  tm_shape(KG_p_africa_poly_shp, projection =  4326) + 
  tm_polygons(col = "#d9d9d9") +
  tm_shape(KG_p_f_arid_NA) + 
  tm_raster(title = "", style = "cat",
            labels = "arid", legend.show = FALSE,
            palette = "#ffffb2") +
  tm_shape(KG_p_f_tropdry_NA) + 
  tm_raster(title = "", style = "cat", legend.show = FALSE,
            labels = "tropical (seasonal)",
            palette = "#762a83") +
  tm_shape(KG_p_f_tropwet_NA) + 
  tm_raster(title = "", style = "cat", legend.show = FALSE,
            labels = "tropical (humid)",
            palette = "#1b7837") +
  tm_shape(KG_p_f_other) + 
  tm_raster(title = "", style = "cat",
            labels = "other", legend.show = FALSE,
            palette = "#fcbba1") +
  tm_shape(KG_p_africa_poly_shp, projection =  4326) + 
  tm_borders(col = "black", lwd = 0.5) +
  tm_shape(Africa_sp, projection = 4326) +
  tm_borders(col = "black", lwd = 0.5) +
  tm_legend(position = c(0.03, 0.05), legend.text.size = 1.4, 
            legend.width = 0.8, legend.title.size = 2, legend.height = 0.8) +
  tm_add_legend(type = "fill", 
                col = c("#d9d9d9"),
                size = 1.5, shape = 0,
                labels = c("Temperate/tropical (seasonal)"),
                title = "Current climate zones") +
  tm_add_legend(type = "fill", 
                col = c("#ffffb2", "#762a83", "#1b7837", "#fcbba1"),
                size = 1.5, shape = 0,
                labels = c("Arid", "Tropical (seasonal)", 
                           "Tropical (humid)", "Other"),
                title = "Predicted changes to") +
  tm_layout(main.title = "a) Predicted changes in climate zones", 
            title.size = 2) 
tmap_save(tm = KG_p_f_NA, 
          filename = "./Figures/AfSIS_14C_Figure4a.jpeg")

##Calculate Area-weighted SOC content
KG_p_area <- as.data.frame(KG_p_africa_grouped) %>% 
  rename(KG = Beck_KG_V1_present_0p0083) %>% 
  filter(KG > 0) %>% 
  mutate(KG_p_group = case_when(
    KG == 1 ~ "Arid",
    KG == 2 ~ "Temperate (seasonal)",
    KG == 3 ~ "Temperate (humid)",
    KG == 4 ~ "Tropical (seasonal)",
    KG == 5 ~ "Tropical (humid)"
  )) %>% 
  group_by(KG_p_group) %>%
  summarise(area = n()) %>%
  mutate(area_all = sum(area),
         rel_area = area / sum(area)) 

KG_p_area$KG_p_group <- factor(KG_p_area$KG_p_group,
                               levels = c("Arid", 
                                          "Temperate (seasonal)",
                                          "Tropical (seasonal)",
                                          "Temperate (humid)",
                                          "Tropical (humid)"))

AfSIS_14C$KG_p_group <- factor(AfSIS_14C$KG_p_group,
                               levels = c("Arid", 
                                          "Temperate (seasonal)",
                                          "Tropical (seasonal)",
                                          "Temperate (humid)",
                                          "Tropical (humid)"))

AfSIS_14C_KG_SOC <- AfSIS_14C %>% 
  group_by(KG_p_group) %>% 
  summarise(mean_se_CORG = mean_se(CORG))

KG_p_area %>%   
  left_join(AfSIS_14C_KG_SOC) %>% 
  ggplot(aes(x = KG_p_group, y = mean_se_CORG$y * rel_area,
             fill = KG_p_group)) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5,
                ymin = 0, ymax = 0.21),
            color = "black", fill = "#d9d9d9") +
  geom_rect(aes(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 0.7),
            color = "black", fill = "#d9d9d9") +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean_se_CORG$ymin * rel_area,
                    ymax = mean_se_CORG$ymax * rel_area),
                width = 0.2) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous("Area-weighted SOC content", expand = c(0,0), 
                     limits = c(0,0.9), breaks = seq(0,0.8,0.2)) +
  scale_x_discrete("", 
                   labels = c("Arid", "Temperate\n(seasonal)", "Tropical\n(seasonal)",
                              "Temperate\n(humid)", "Tropical\n(humid)")) +
  scale_fill_manual("", values = c("#ffffb2",
                                   "#c2a5cf", "#762a83", 
                                   "#a6dba0", "#1b7837")) +
  guides(fill = "none") +
  labs(title = "b) Mean SOC content of current climate zones")
ggsave("./Figures/AfSIS_14C_Figure4b.jpeg", width = 5.6, height = 2.5) 
