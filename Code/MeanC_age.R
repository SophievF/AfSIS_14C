## Calculate mean C age from 14C values from the AfSIS dataset ##
## Sophie von Fromm ##
## 2022-04-16 ##

## Calculate mean C age based on an one-pool model from 14C data
# Code adopted from (Khomo et al. 2017, SI p. 19) and with input from Jeff Beem-Miller

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

# Define atmospheric 14C curve for the different regions/records
NHZone3 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$NHZone3, 
                          time.scale = "AD")
SHZone12 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$SHZone12, 
                           time.scale = "AD")
SHZone3 <- bind.C14curves(prebomb = IntCal13, postbomb = Hua2013$SHZone3, 
                          time.scale = "AD")

# Define model parameters for calculation 
ModelParameter_TT <- data.frame(
  TT = seq_len(6500), ## Turnover time to search for matching 14C value
  In = 100) %>%  #Arbitrary value since system is in steady state
  mutate(k = 1/TT, #Decomposition rate: 1/turnover time in 1/yr
         F0 = k/(k+(1/8267)), #Fraction modern in the soil under steady state conditions
         DeltaF0 = (F0-1) * 1000, #Calculate Delta14C from F0
         C0 = In * TT) #Carbon stock in the soil Cstock = Input/k

## Function to run model
# k: column 3, C0: column 6, F0 (in Delta14C): column 5, In: column 2
opm_fun <- function(zone){
  ModelParameter_TT %>% 
    pmap(~OnepModel14(t = seq(1901, 2009, by = 0.5), 
                      k = ..3, C0 = ..6,
                      F0 = ..5, In = ..2, 
                      inputFc = zone))
}

## Run model for each zone
# NHZone3
opm_NHZ3 <- opm_fun(zone = NHZone3)

#SHZone12
opm_SHZ12 <- opm_fun(zone = SHZone12)

#SHZone3
opm_SHZ3 <- opm_fun(zone = SHZone3)

## Function to extract 14C for year 2009 (last entry in model list = 217)
# takes quite some time to run model
d14C_extract <- function(model){
  model %>% 
    map(getF14) %>% 
    map(217)
}

# output: list for each turnover time (1-6500yrs) for the year 2009

#NHZone3
opm_d14C_2009_NHZ3 <- d14C_extract(model = opm_NHZ3)

#SHZone12
opm_d14C_2009_SHZ12 <- d14C_extract(model = opm_SHZ12)

#SHZone3
opm_d14C_2009_SHZ3 <- d14C_extract(model = opm_SHZ3)

## Function to convert list into a data frame to check for highest 14C value
# Due to "bomb14C" there are two possible turnover times for high 14C values
opm_df <- function(model_results){
  data.frame(
    Delta14C_mod = unlist(model_results),
    TurnoverTime = ModelParameter_TT$TT
  )
}

#NHZone3
opm_NHZ3_df <- opm_df(model_results = opm_d14C_2009_NHZ3)

#SHZone12
opm_SHZ12_df <- opm_df(model_results = opm_d14C_2009_SHZ12)

#SHZone3
opm_SHZ3_df <- opm_df(model_results = opm_d14C_2009_SHZ3)

## Position of highest 14C value (bomb peak)
# only consider higher values to make sure that always the longer TT get's selected
# older turnover times are more realistic for bulk 14C measurements
which.max(opm_NHZ3_df$Delta14C_mod) 
which.max(opm_SHZ12_df$Delta14C_mod) 
which.max(opm_SHZ3_df$Delta14C_mod) 

#Function to extract observed/measured 14C values for each zone
d14C_obs <- function(zone){
  df_14c_zones %>% 
    filter(Zones == zone) %>% 
    .$Delta14C
}

#NHZone3
d14C_obs_NHZ3 <- d14C_obs(zone = "NHZone3")

#SHZone12
d14C_obs_SHZ12 <- d14C_obs(zone = "SHZone12")

#SHZone3
d14C_obs_SHZ3 <- d14C_obs(zone = "SHZone3")

# Function to create empty list to find closest match between observed and modeled 14C
empty_list <- function(data){
  vector(mode = "list", length = length(data))
}

#NHZone3
d14C_obs_NHZ3_l <- empty_list(data = d14C_obs_NHZ3)

#SHZone12
d14C_obs_SHZ12_l <- empty_list(data = d14C_obs_SHZ12)

#SHZone3
d14C_obs_SHZ3_l <- empty_list(data = d14C_obs_SHZ3)

## Find closest match between observed an modeled 14C

#NHZone3
d14C_obs_NHZ3_l <- map(seq_along(d14C_obs_NHZ3), function(i){
  d14C_obs_NHZ3_l[i] <- which.min(abs(unlist(opm_d14C_2009_NHZ3[27:6500]) - d14C_obs_NHZ3[i]))
})

#SHZone12
d14C_obs_SHZ12_l <- map(seq_along(d14C_obs_SHZ12), function(i){
  d14C_obs_SHZ12_l[i] <- which.min(abs(unlist(opm_d14C_2009_SHZ12[26:6500]) - d14C_obs_SHZ12[i]))
})

#SHZone3
d14C_obs_SHZ3_l <- map(seq_along(d14C_obs_SHZ3), function(i){
  d14C_obs_SHZ3_l[i] <- which.min(abs(unlist(opm_d14C_2009_SHZ3[26:6500]) - d14C_obs_SHZ3[i]))
})


## Extract and merge results

#NHZone3
TT_NHZ3_df <- data.frame(Delta14C = unlist(d14C_obs_NHZ3), 
                         TurnoverTime = unlist(d14C_obs_NHZ3_l),
                         Zones = "NHZone3")

#SHZone12
TT_SHZ12_df <- data.frame(Delta14C = unlist(d14C_obs_SHZ12), 
                          TurnoverTime = unlist(d14C_obs_SHZ12_l),
                          Zones = "SHZone12")

#SHZone3
TT_SHZ3_df <- data.frame(Delta14C = unlist(d14C_obs_SHZ3), 
                         TurnoverTime = unlist(d14C_obs_SHZ3_l),
                         Zones = "SHZone3")



TT_SHZ3_df %>% 
  ggplot(aes(x = Delta14C, y = TurnoverTime)) +
  geom_point()

# Merge all three data frames and add back to original data
df_14C_TT <- rbind(TT_NHZ3_df, TT_SHZ12_df, TT_SHZ3_df) %>% 
  right_join(df_14c_zones, by = c("Delta14C", "Zones")) %>%
  #remove duplicate rows that are created due to samples that have the same 14C value within one zone
  unique() %>% 
  tibble()

# Plot results
df_14C_TT %>% 
  ggplot(aes(x = Delta14C, y = TurnoverTime, color = Zones)) +
  geom_point()

write.csv(df_14C_TT, "./Data/df_14C_TT.csv", row.names = FALSE)
