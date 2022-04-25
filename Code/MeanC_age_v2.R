## Sophie von Fromm
## 2022-April-25

## Calculate mean C age based on an one-pool model from 14C data

# Code adopted from (Khomo et al. 2017, SI p. 19) and with input from Jeff Beem-Miller
# Idea for code from Sierra et al. 2014 and radiocarbon book (chapter 3)

## Load SoilR
installed_packages <- "SoilR" %in% rownames(installed.packages())

if (any(installed_packages == FALSE)){
  devtools::install_github("MPIBGC-TEE/SoilR-exp/pkg",
                           upgrade = FALSE)
}

library(SoilR)
library(tidyverse)

# Define zones for 14C atmosphere curves (zones defined on latitudes)
# Load 14C data and long/lat data
source("./Code/Load14C_Data.R")

tbl_14c <- read_csv_14C_fun(path = "./Data/Radiocarbon")

tbl_longlat <- read_csv("./Data/AfSIS_LongLat.csv")

df_14c <- tbl_longlat %>% 
  full_join(tbl_14c)

zones_14C <- function(x){
  x %>% 
    mutate(Zones = case_when(
      Latitude > 0 ~ "NHZone3",
      Latitude < 0 & Latitude > -15 ~ "SHZone12",
      Latitude < -15 ~ "SHZone3"
    ))
}

df_14c_zones <- zones_14C(df_14c)

#Check Zones
df_14c_zones %>% 
  group_by(Zones) %>% 
  count()

