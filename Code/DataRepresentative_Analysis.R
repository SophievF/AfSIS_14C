## Sophie von Fromm
## 2022-November-07

library(tidyverse)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(ISRaD)
library(sf)
library(raster)
library(grid)
library(tmap)
library(ncdf4)
library(ISRaD)

#This script contains all the code to reproduce all data analysis and figures 
#in the manuscript von Fromm et al. (2023)

#Load dataset
AfSIS_14C <- read_csv("./Data/AfSIS_data_all.csv") %>% 
  dplyr::rename(Cultivation = PlotCultMgd)

AfSIS_14C$KG_p_group <- factor(AfSIS_14C$KG_p_group,
                               levels = c("Arid", 
                                          "Temperate (seasonal)",
                                          "Tropical (seasonal)",
                                          "Temperate (humid)",
                                          "Tropical (humid)")) 

AfSIS_14C$Depth <- factor(AfSIS_14C$Depth, 
                          levels = c("Topsoil", "Subsoil"))

AfSIS_14C$Erosion <- factor(AfSIS_14C$Erosion)

AfSIS_14C$Cultivation <- factor(AfSIS_14C$Cultivation)

##Preparation for figures
#Create own scale for mean C age
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#Create own theme
theme_own <- theme(axis.text = element_text(color = "black"),
                   plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   strip.background = element_rect(fill = NA),
                   strip.text = element_text(face = "bold", size = 14),
                   axis.title = element_text(face = "bold"))

names(AfSIS_14C)

##### Sampling density #####

ISRaD_dir <- "C:/Users/sfromm/Documents/GitHub/ISRaD/ISRaD_data_files"
ISRaD_extra <- ISRaD.getdata(directory = ISRaD_dir,
                             dataset = "full", extra = TRUE,
                             force_download = TRUE)

ISRaD_lyr <- ISRaD.flatten(ISRaD_extra, "layer") %>% 
  drop_na(lyr_14c) %>% 
  # filter(is.na(pro_land_cover) | pro_land_cover != "wetland") %>% 
  # filter(is.na(pro_usda_soil_order) | pro_usda_soil_order != "Histosols") %>% 
  tibble()

ISRaD_lyr %>% 
  count(entry_name)

ISRaD_lyr %>% 
  count(pro_land_cover) %>% 
  mutate(rel_dis = n/sum(n)*100)

ISRaD_lyr %>%  
  mutate(usda_soil = case_when(
    is.na(pro_usda_soil_order) ~ pro_soilOrder_0.5_deg_USDA.1,
    TRUE ~ pro_usda_soil_order
  )) %>% 
  count(usda_soil) %>% 
  mutate(rel_dis = n/sum(n)*100)

data("World")

world_map <- World %>% 
  filter(continent != "Seven seas (open ocean)")

world_map_w_SSA <- world_map %>% 
  filter(name != "Algeria",
         name != "Egypt",
         name != "Libya",
         name != "W. Sahara",
         name != "Tunisia",
         name != "Morocco",
         name != "Djibouti",
         name != "Sudan")

world_map$continent <- replace(world_map$continent,
                               which(world_map$name == "Russia"),
                               "Asia")
world_map_w_SSA$continent <- replace(world_map_w_SSA$continent,
                                      which(world_map_w_SSA$name == "Russia"),
                                      "Asia")

area_world_SSA <- world_map_w_SSA %>% 
  as.data.frame() %>% 
  group_by(continent) %>% 
  summarise(area_sum = sum(area))

st_crs(world_map) <- 4326
st_crs(world_map_w_SSA) <- 4326

ISRaD_sf <- sf::st_as_sf(ISRaD_lyr, coords = c("pro_long", "pro_lat"), crs = 4326)

afsis_sf <- sf::st_as_sf(AfSIS_14C, coords = c("Longitude", "Latitude"), crs = 4326)

map_sf <- sf::st_join(ISRaD_sf, world_map_w_SSA)

map_data <- st_set_geometry(map_sf, NULL)

map_all <- map_data %>% 
  full_join(ISRaD_lyr) %>% 
  tibble()

world_map_df <- data.frame(as_Spatial(world_map))

# Map global and AfSIS sampling distribution
map14c <- ggplot() +
  geom_sf(data = world_map, fill = NA) +
  geom_sf(data = world_map_w_SSA, aes(fill = continent),
          alpha = 0.5) +
  geom_jitter(data = map_all, aes(x = pro_long, y = pro_lat),
              size = 2, width = 0.5, height = 0.5, shape = 1) +
  geom_jitter(data = AfSIS_14C, aes(x = Longitude, y = Latitude),
              color = "red", size = 2, width = 0.5, height = 0.5, shape = 2) +
  geom_rect(aes(ymin = -90, ymax = -60, xmin = -180, xmax = 0),
            fill = "white", alpha = 0.9) +
  theme_bw(base_size = 14) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), expand = c(0,0),
                     breaks = c(-100,0,100), limits = c(-180,180)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), expand = c(0,0),
                     breaks = c(-50,0,50), limits = c(-90,90)) +
  # scale_fill_brewer(palette = "YlGnBu") +
  scale_fill_manual(values = c("#ffffcc", "lightgrey", "lightgrey", "lightgrey",
                               "lightgrey", "lightgrey", "lightgrey"))

#Gap-fill missing continent
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "Antarctica" &
                                     is.na(map_all$continent)),
                             "Antarctica")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "Canada" &
                                     is.na(map_all$continent)),
                             "North America")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "France" &
                                     is.na(map_all$continent)),
                             "Europe")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "Italy" &
                                     is.na(map_all$continent)),
                             "Europe")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "Malaysia" &
                                     is.na(map_all$continent)),
                             "Asia")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "Panama" &
                                     is.na(map_all$continent)),
                             "North America")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "Russia" &
                                     is.na(map_all$continent)),
                             "Asia")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_country == "South Africa" &
                                     is.na(map_all$continent)),
                             "Africa")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_name == "Pu'u Eke 1" &
                                     is.na(map_all$continent)),
                             "North America")
map_all$continent <- replace(map_all$continent,
                             which(map_all$entry_name == "Mergelov_2020" &
                                     is.na(map_all$continent)),
                             "Antarctica")
map_all$continent <- replace(map_all$continent,
                             which(map_all$entry_name == "Zazovskaya_2017" &
                                     is.na(map_all$continent)),
                             "Antarctica")
map_all$continent <- replace(map_all$continent,
                             which(map_all$site_name == "GORG" &
                                     is.na(map_all$continent)),
                             "South America")
map_all$continent <- replace(map_all$continent,
                             which(map_all$site_name == "Molino Meloni" &
                                     is.na(map_all$continent)),
                             "Europe")
map_all$continent <- replace(map_all$continent,
                             which(map_all$site_name == "Molino Meloni" &
                                     is.na(map_all$continent)),
                             "Europe")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_name == "Nurallao:39.78,8.38" &
                                     is.na(map_all$continent)),
                             "Europe")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_name == "Grab 89-14" &
                                     is.na(map_all$continent)),
                             "Asia")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_name == "Grab 89-16" &
                                     is.na(map_all$continent)),
                             "Asia")
map_all$continent <- replace(map_all$continent,
                             which(map_all$entry_name == "Staub_2003" &
                                     is.na(map_all$continent)),
                             "Asia")
map_all$continent <- replace(map_all$continent,
                             which(map_all$pro_name == "TCC2" &
                                     is.na(map_all$continent)),
                             "South America")

map_all %>% 
  summarise(n_study = n_distinct(entry_name),
            n = n())

map_all %>% 
  count(continent)

map_all %>% 
  filter(continent == "Africa") %>% 
  count(entry_name)

map_all %>% 
  summarise(n = n()/n_distinct(entry_name))

# Calculate sampling density
dens_14c <- map_all %>% 
  group_by(continent) %>% 
  summarise(n_14c = n()) %>% 
  mutate(n_14c_afsis = case_when(
    continent == "Africa" ~ n_14c + 514)) %>% 
  right_join(area_world_SSA) %>% 
  mutate(dens_all = n_14c/area_sum*1000000,
         dens_ssa_new = n_14c_afsis/area_sum*1000000)

# plot sampling density
dens_bar <- dens_14c %>% 
  ggplot() +
  geom_bar(aes(x = reorder(continent, -as.numeric(dens_all)), 
               y = as.numeric(dens_ssa_new), fill = "AfSIS"),
           stat = "identity", color = "black") +
  geom_bar(aes(x = reorder(continent, -as.numeric(dens_all)), 
               y = as.numeric(dens_all), fill = "ISRaD"), 
           stat = "identity", color = "black") +
  theme_classic(base_size = 12) +
  theme_own +
  theme(legend.position = c(0.42, 0.8),
        panel.background = element_blank(),
        legend.background = element_blank()) +
  scale_x_discrete("", labels = c("Europe", "N. America", "Oceania",
                                  "S. America", "Asia", "SSA", "Antartica")) +
  scale_y_continuous(expression("Samples / km² * 10"^-6), expand = c(0,0),
                     limits = c(0,200)) +
  scale_fill_manual("Dataset", values = c("red", "lightgrey")) 

# plot both together
map14c_densbar <- map14c +
  annotation_custom(ggplotGrob(dens_bar), xmin = -180, xmax = 0,
                    ymin = -93, ymax = 5)
ggsave("./Figures/AfSIS_14C_Figure1.jpeg", width = 12, height = 7)

##### Data distribution/representative #####

##Reference data can be downloaded from: https://data.worldagroforestry.org/dataset.xhtml?persistentId=doi:10.34725/DVN/66BFOB
afsis_ref <- read_csv("H:/PhD/AfSIS/R_Data/AfSIS_RefData_Roth_allSites_GlobalData.csv") %>% 
  mutate(Alox = Am_Ox_Al/10000,
         Feox = Am_Ox_Fe/10000) 

AfSIS_GPP <- read_csv("./Data/AfSIS_GPP.csv") %>% 
  dplyr::rename(SSN = SiteID) %>% 
  #Coordinates have one digit less than original data
  dplyr::select(SSN, GPP)

afsis_14c <- read_csv("./Data/df_14C_TT.csv") %>% 
  dplyr::select(SSN, Delta14C)

afsis_all <- afsis_ref %>% 
  left_join(AfSIS_GPP) %>% 
  left_join(afsis_14c) %>% 
  mutate(data_14c = case_when(
    is.na(Delta14C) ~ FALSE,
    TRUE ~ TRUE
  )) %>% 
  mutate(data_ref = TRUE) %>% 
  dplyr::select(SSN, Longitude, Latitude, Region, Country, Site, Cluster, Profile,
                Depth, Clay_8um, CORG, pH, Alox, Feox, MAT, MAP, AridityIndex, GPP,
                data_14c, data_ref) 

afsis_all %>% 
  dplyr::count(data_ref)

## plot data distribution for AfSIS 14C and AfSIS reference
# https://aosmith.rbind.io/2018/08/20/automating-exploratory-plots/

feat <- names(afsis_all)[10:18]
feat <- set_names(feat)

fun_plot <- function(x){
  ggplot(mapping = aes(x = .data[[x]])) +
    geom_density(data = afsis_all, aes(color = "black"), trim = "TRUE") +
    geom_density(data = afsis_all %>% 
                   filter(data_14c == TRUE), aes(color = data_14c), trim = "TRUE") +
    theme_classic(base_size = 18) +
    theme_own +
    scale_y_continuous("", labels = label_comma(accuracy = 0.0001), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_manual("Dataset", labels = c("AfSIS_ref", "AfSIS_14c"),
                       values = c("black", "red"))
}

all_plots <- map(feat, ~fun_plot(.x))

ggarrange(plotlist = all_plots, common.legend = TRUE)
ggsave("./Figures/AfSIS_14C_FigureA4.jpeg", width = 12, height = 6)


### Distribution of selected parameters for AfSIS 14C, ref and SSA

## Climate data

# MAT
MAT_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_01.tif"
MAT_global <- raster(MAT_dir)

# MAP
MAP_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_12.tif"
MAP_global <- raster(MAP_dir)

# PET
PET_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/global-et0_annual.tif/et0_yr/et0_yr.tif"
PET_global <- raster(PET_dir)  

## GPP (from Fluxcom)
# GPP from FLUXCOM
read_nc_fun <- function(path, pattern = "*.nc"){
  list.files(path, pattern, full.names = TRUE) %>% 
    map(~brick(.))
}

GPP_global <- read_nc_fun(path = "D:/Sophie/PhD/AfSIS_GlobalData/Fluxcom_GPP")

#calculate yearly means for each year
Global_calc <- map(GPP_global, ~calc(., fun = mean))

#Calculate long-term mean
GPP_mean <- (Global_calc[[1]]+Global_calc[[2]]+Global_calc[[3]]+
               Global_calc[[4]]+Global_calc[[5]]+Global_calc[[6]]+
               Global_calc[[7]]+Global_calc[[8]]+Global_calc[[9]]+
               Global_calc[[10]]+Global_calc[[11]]+Global_calc[[12]])/12

## Soil properties
# pH (Hengl et al 2021; 30m resolution)
pH_dir_top <- "D:/Sophie/PhD/AfSIS_GlobalData/Africa_SoilMaps_30m_2022/sol_ph_h2o_m_30m_0..20cm_2001..2017_v0.13_wgs84.tif"
pH_dir_bot <- "D:/Sophie/PhD/AfSIS_GlobalData/Africa_SoilMaps_30m_2022/sol_ph_h2o_m_30m_20..50cm_2001..2017_v0.13_wgs84.tif"
pH_all <- stack(pH_dir_top, pH_dir_bot)

# Boundaries for Africa
data("World")
africa_map <- World %>%
  st_transform(4326) %>% 
  dplyr::filter(continent == "Africa")

# Sub-Saharan Africa (based on definition of the United Nations geoscheme for Africa)
ssa_sp <- africa_map %>% 
  filter(name != "Algeria",
         name != "Egypt",
         name != "Libya",
         name != "W. Sahara",
         name != "Tunisia",
         name != "Morocco",
         name != "Djibouti",
         name != "Sudan") %>% 
  as("Spatial")

# Crop climate data to boundaries for SSA
MAT_crop <- raster::crop(MAT_global, ssa_sp)
MAT_africa <- raster::mask(MAT_crop, ssa_sp)

MAP_crop <- raster::crop(MAP_global, ssa_sp)
MAP_africa <- raster::mask(MAP_crop, ssa_sp)

PET_crop <- raster::crop(PET_global, ssa_sp)
PET_africa <- raster::mask(PET_crop, ssa_sp)  
  
AI_africa <- PET_africa/MAP_africa

# Crop GPP to boundaries for SSA
GPP_crop <- raster::crop(GPP_mean, ssa_sp)
GPP_africa <- raster::mask(GPP_crop, ssa_sp)

## Randomly selected 2,002 sampling points across SSA
set.seed(42)
sf_use_s2(FALSE)

rdm_p <- sampleRandom(MAP_africa, size = 2002, sp = TRUE)

plot(rdm_p)

# MAP
MAP_rdm <- raster::extract(MAP_africa, rdm_p, df = TRUE) %>% 
  rename(MAP = wc2.0_bio_30s_12)

# MAT
MAT_rdm <- raster::extract(MAT_africa, rdm_p, df = TRUE) %>% 
  rename(MAT = wc2.0_bio_30s_01)

#GPP
GPP_rdm <- raster::extract(GPP_africa, rdm_p, df = TRUE) %>% 
  rename(GPP = layer)

# pH
pH_rdm <- raster::extract(pH_all, rdm_p, df = TRUE) %>% 
  rename(d_0_20 = sol_ph_h2o_m_30m_0..20cm_2001..2017_v0.13_wgs84,
         d_20_50 = sol_ph_h2o_m_30m_20..50cm_2001..2017_v0.13_wgs84) %>% 
  pivot_longer(!ID, names_to = "depth", values_to = "pH")

rdm_data <- MAT_rdm %>% 
  full_join(MAP_rdm) %>% 
  full_join(pH_rdm) %>% 
  full_join(GPP_rdm) %>% 
  tibble() %>% 
  mutate(data = "SSA") %>% 
  mutate(pH = pH/10)

skimr::skim_without_charts(rdm_data)

# Violin plots: SSA and AfSIS 14C
p1 <- afsis_all %>%
  filter(data_14c == TRUE) %>% 
  dplyr::select(GPP, MAP, MAT, pH, data_14c) %>% 
  full_join(rdm_data %>% 
              mutate(data_14c = FALSE) %>% 
              dplyr::select(GPP, MAT, pH, MAP, data_14c)) %>% 
  ggplot() + 
  geom_violin(aes(y = GPP*365/1000, x = data_14c, fill = data_14c),
              draw_quantiles = c(0.25,0.75)) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS")) +
  scale_y_continuous("GPP [kgC/m²yr]") +
  scale_fill_manual(values = c("#ffffcc", "#FF7F7F"))

p2 <- afsis_all %>%
  filter(data_14c == TRUE) %>% 
  dplyr::select(GPP, MAP, MAT, pH, data_14c) %>% 
  full_join(rdm_data %>% 
              mutate(data_14c = FALSE) %>% 
              dplyr::select(GPP, MAT, pH, MAP, data_14c)) %>% 
  ggplot() + 
  geom_violin(aes(y = MAP, x = data_14c, fill = data_14c),
              draw_quantiles = c(0.25,0.75)) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS")) +
  scale_y_continuous("MAP [mm]", expand = c(0,0), limits = c(0,3600)) +
  scale_fill_manual(values = c("#ffffcc", "#FF7F7F")) 

p3 <- afsis_all %>%
  filter(data_14c == TRUE) %>% 
  dplyr::select(GPP, MAP, MAT, pH, data_14c) %>% 
  full_join(rdm_data %>% 
              mutate(data_14c = FALSE) %>% 
              dplyr::select(GPP, MAT, pH, MAP, data_14c)) %>% 
  ggplot() + 
  geom_violin(aes(y = pH, x = data_14c, fill = data_14c),
              draw_quantiles = c(0.25,0.75)) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS")) +
  scale_y_continuous("pH", expand = c(0,0), limits = c(4,10)) +
  scale_fill_manual(values = c("#ffffcc", "#FF7F7F")) 

p4 <- afsis_all %>%
  filter(data_14c == TRUE) %>% 
  dplyr::select(GPP, MAP, MAT, pH, data_14c) %>% 
  full_join(rdm_data %>% 
              mutate(data_14c = FALSE) %>% 
              dplyr::select(GPP, MAT, pH, MAP, data_14c)) %>% 
  ggplot() + 
  geom_violin(aes(y = MAT, x = data_14c, fill = data_14c),
              draw_quantiles = c(0.25,0.75)) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS")) +
  scale_y_continuous("MAT [°C]", expand = c(0,0), limits = c(5,33)) +
  scale_fill_manual(values = c("#ffffcc", "#FF7F7F")) 

ggarrange(p1,p2,p3,p4)
ggsave("./Figures/AfSIS_14C_Figure1c.jpeg", 
       width = 6, height = 6)


# Data distribution SSA, AfSIS ref, AfSIS 14C
p_gpp <- ggplot() +
  geom_violin(data = afsis_all, mapping = aes(y = GPP, x = data_14c, fill = data_14c), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = afsis_all, mapping = aes(y = GPP, x = data_14c), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 1/20) +
  geom_violin(data = GPP_rdm,  mapping = aes(y = GPP, x = "blue", fill = "blue"), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = GPP_rdm,  mapping = aes(y = GPP, x = "blue"), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 1/20) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS ref","AfSIS 14C"), expand = c(0,0)) +
  scale_y_continuous("GPP [gC/m²d]") +
  scale_fill_brewer(palette = "Paired") 

p_MAP <- ggplot() +
  geom_violin(data = afsis_all, mapping = aes(y = MAP, x = data_14c, fill = data_14c), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = afsis_all, mapping = aes(y = MAP, x = data_14c), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 20) +
  geom_violin(data = MAP_rdm,  mapping = aes(y = MAP, x = "blue", fill = "blue"), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = MAP_rdm,  mapping = aes(y = MAP, x = "blue"), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 20) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS ref","AfSIS 14C"), expand = c(0,0)) +
  scale_y_continuous("MAP [mm]") +
  scale_fill_brewer(palette = "Paired")  

p_MAT <- ggplot() +
  geom_violin(data = afsis_all, mapping = aes(y = MAT, x = data_14c, fill = data_14c), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = afsis_all, mapping = aes(y = MAT, x = data_14c), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 1/10) +
  geom_violin(data = MAT_rdm,  mapping = aes(y = MAT, x = "blue", fill = "blue"), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = MAT_rdm,  mapping = aes(y = MAT, x = "blue"), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 1/10) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS ref","AfSIS 14C"), expand = c(0,0)) +
  scale_y_continuous("MAT [°C]") +
  scale_fill_brewer(palette = "Paired") 

p_pH <- ggplot() +
  geom_violin(data = afsis_all, mapping = aes(y = pH, x = data_14c, fill = data_14c), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = afsis_all, mapping = aes(y = pH, x = data_14c), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 1/25) +
  geom_violin(data = pH_rdm,  mapping = aes(y = pH/10, x = "blue", fill = "blue"), 
              width = 0.9, draw_quantiles = c(0.25,0.75)) +
  geom_dotplot(data = pH_rdm,  mapping = aes(y = pH/10, x = "blue"), binaxis = "y",
               stackdir = "center", dotsize = 0.3, binwidth = 1/25) +
  theme_classic(base_size = 16) +
  theme_own +
  theme(legend.position = "none") +
  scale_x_discrete("",labels = c("SSA", "AfSIS ref","AfSIS 14C"), expand = c(0,0)) +
  scale_y_continuous("pH") +
  scale_fill_brewer(palette = "Paired") 

ggarrange(p_gpp, p_MAP, p_MAT, p_pH)
ggsave("./Figures/AfSIS_14C_FigureA2.jpeg", width = 12, height = 6)

             