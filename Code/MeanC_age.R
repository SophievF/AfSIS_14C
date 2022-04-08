## Sophie von Fromm
## 2022-April-08

## Calculate mean C age based on an one-pool model from 14C data

# Code adopted from Sue (Khomo et al. 2017, SI p. 19) and Jeff
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
      Latitude < 0 & Latitude -15 ~ "SHZone12",
      Latitude < -15 ~ "SHZone3"
    ))
}

df_14c_zones <- zones_14C(df_14c)

# Define 14C atmosphere curves for different zones
curve_list <- list(
  NHZone3 = NHZone3 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$NHZone3, 
                                      time.scale = "AD"), 
  SHZone12 = SHZone12 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$SHZone12, 
                                        time.scale = "AD"),
  SHZone3 = SHZone3 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$SHZone3, 
                                      time.scale = "AD"))


# Define model parameters
model_par <- data.frame(
  TT = seq_len(6500), # possible turnover times to search for matching 14C values
  In = 100) %>% # arbitrary C input value since system is assumed to be in steady state
  mutate(k = 1/TT, # decomposition rate: 1/turn over time in 1/yr
         F0 = k/(k+(1/8267)), # fraction modern in the soil under steady state conditions
         DeltaF0 = (F0-1)*100, # calculate delta14C from fraction modern
         C0 = In * TT) # carbon stock in th soil; Cstock = Input/k

# Run one pool model for each possible turnover and for the three different zones
# k - column 3, C0 - column 6, DeltaF0 - column 5, In -column 2
feat <- names(curve_list)

#Working
one_pool_model <- function(x){
  model_par %>% 
  pmap(~OnepModel14(t = seq(1901,2009, by = 0.5),
                    k =  ..3, C0 = ..6,
                    F0_Delta14C = ..5, 
                    In =..2, inputFc = curve_list[[x]]))
}

#not really tested; change assignment: use variable names instead of column number?!
one_pool_model <- function(x){
  pmap(model_par, ~OnepModel14(t = seq(1901,2009, by = 0.5),
                               k =  ..3, C0 = ..6,
                               F0_Delta14C = ..5, 
                               In =..2, inputFc = curve_list[[x]]))
}

#not tested
one_pool_model <- function(x){
  pmap(model_par, ~getF14(OnepModel14(t = seq(1901,2009, by = 0.5),
                                      k =  ..3, C0 = ..6,
                                      F0_Delta14C = ..5, 
                                      In =..2, inputFc = curve_list[[x]])))
}


# Try not to store opm because it is so large

opm <- map(feat, ~one_pool_model(.x))


# Extract 14C fraction for each zone
#not working: opm_F14 <- map(opm, getF14)

#not ideal: need to do everystep three times...
opm_NHZ3 <- opm[[1]] %>% 
  map(getF14) 

# Extract values for year 2009 (last entry in list)
opm_NHZ3_2009 <- opm_NHZ3 %>% 
  map(217) 
  
df_NHZ3_2009 <- data.frame(
    Delta14C_TT = unlist(opm_NHZ3_2009),
    TT = model_par$TT)



opm_SHZ12 <- opm[[2]] %>% 
  map(getF14)
opm_SHZ3 <- opm[[3]] %>% 
  map(getF14)


