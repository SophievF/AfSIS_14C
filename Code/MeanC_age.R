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
      Latitude < 0 & Latitude > -15 ~ "SHZone12",
      Latitude < -15 ~ "SHZone3"
    ))
}

df_14c_zones <- zones_14C(df_14c)

#Check Zones
df_14c_zones %>% 
  group_by(Zones) %>% 
  count()


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
  #Use smaller number for testing
  TT = seq_len(6500), # possible turnover times to search for matching 14C values
  In = 100) %>% # arbitrary C input value since system is assumed to be in steady state
  mutate(k = 1/TT, # decomposition rate: 1/turn over time in 1/yr
         F0 = k/(k+(1/8267)), # fraction modern in the soil under steady state conditions
         DeltaF0 = (F0-1)*100, # calculate delta14C from fraction modern
         C0 = In * TT) # carbon stock in the soil; Cstock = Input/k

# Run one pool model for each possible turnover and for the three different zones
# k - column 3, C0 - column 6, DeltaF0 - column 5, In -column 2


#Function to run model for each atmospheric zone (n = 3)
feat <- names(curve_list)

one_pool_model <- function(x){
  pmap(model_par, ~OnepModel14(t = seq(1901,2009, by = 0.5),
                               k =  ..3, C0 = ..6,
                               F0_Delta14C = ..5, 
                               In =..2, inputFc = curve_list[[x]]))
}

# function to apply model and extract values for target year (2009) by zone
opm_fun <- function(x){
  opm <- map(feat, ~one_pool_model(.x))
  names(opm) <- feat
  opm_2009 <- opm[[x]] %>% 
    map(getF14) %>% 
    map(217)
}

# create list that contains model results (14C) for each atmospheric zone
# takes a lot time to create list > 20min
opm_2009_list <- map(feat, ~opm_fun(.x))

# function to create list that contains a data.frame for each zone with 14C and TT
opm_list_fun <- function(x){
  opm_df <- list(data.frame(
    Delta14C_TT = unlist(x),
    TT = model_par$TT
  ))
}

opm_list <- map(opm_2009_list, ~opm_list_fun(.x))
names(opm_list) <- feat


#Find index of maximum value in each list
max_ind <- map(opm_list, function(x){
  which.max(x[[1]]$Delta14C_TT)
})

#NHZ3
max_ind[[1]]


##STOPPED HERE

#Make figures that compare measured and modelled 14C values
  

## Extract TT where difference between modelled and measured 14C is smallest
#Extract measured 14C values for NHZ3
d14c_obs_NHZ3 <- df_test %>% 
  filter(Zones == "NHZone3") %>% .$Delta14C  # observed 14C values for NHZ3

# list of indices with closest match to d14c_obs
d14c_obs_NHZ3_ix <- vector(mode = "list", length = length(d14c_obs_NHZ3))
d14c_obs_NHZ3_ix <- lapply(seq_along(d14c_obs_NHZ3), function(i) {
  d14c_obs_NHZ3_ix[i] <- which.min(abs(unlist(opm_NH3_f14_2009[27:2000]) - d14c_obs_NHZ3[i]))
}) 

TT <- seq(1,2000)
TurnoverTime <- vector(mode = "list", length = length(d14c_obs_NHZ3_ix))
TurnoverTime <- lapply(seq_along(d14c_obs_NHZ3_ix), function(i) {
  TurnoverTime[i] <- TT[d14c_obs_NHZ3_ix[[i]]]
}) 

# compare output
AfSIS_14C_TurnoverTime_NHZ3 <- data.frame(Delta14C = unlist(d14c_obs_NHZ3), 
                                          TurnoverTime = unlist(TurnoverTime),
                                          Zones = "NHZone3")

view(AfSIS_14C_TurnoverTime_NHZ3)




