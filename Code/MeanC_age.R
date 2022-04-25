## Sophie von Fromm
## 2022-April-08


###CODE IS NOT WORKING####

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

# apply both functions at the same time
# takes a lot time to create list > 20min
opm_list <- map(
  map(feat, ~opm_fun(.x)), ~opm_list_fun(.x)
)
names(opm_list) <- feat

#Find index of maximum value in each list
max_ind <- map(opm_list, function(x){
  which.max(x[[1]]$Delta14C_TT)
})


## Extract TT where difference between modeled and measured 14C is smallest
# Function to extract observed 14C values for each zone
d14C_measured_fun <- function(x){
  df_14c_zones %>% 
    filter(Zones == x) %>% 
   .$Delta14C
}

d14C_obs <- map(feat, ~d14C_measured_fun(.x))

# list of indices with closest match to d14c_obs
empty_list_fun <- function(x){
  vector(mode = "list", length = length(x))
}

d14C_obs_NHZ3_ix <- empty_list_fun(x = d14C_obs[[1]])
d14C_obs_SHZ12_ix <- empty_list_fun(x = d14C_obs[[2]])
d14C_obs_SHZ3_ix <- empty_list_fun(x = d14C_obs[[3]])


# function to extract TT where observed and measured values are closest
d14C_match_fun <- function(i, y, x, z) {
  x[i] <- which.min(abs(y - z[i]))
}

# apply function for each zone and store results in data.frame
d14C_obs_NHZ3_ix <- map(seq_along(seq_along(d14C_obs[[1]])), 
                        ~d14C_match_fun(i = .x, x = d14C_obs_NHZ3_ix, 
                                        y = unlist(opm_list$NHZone3)[max_ind[[1]]:6500], 
                                        z = d14C_obs[[1]]))


df_TT_NHZ3 <- data.frame(TurnoverTime = unlist(d14C_obs_NHZ3_ix),
                         Delta14C = unlist(d14C_obs[[1]]),
                         Zones = feat[1])

d14C_obs_SHZ12_ix <- map(seq_along(seq_along(d14C_obs[[2]])), 
                         ~d14C_match_fun(i = .x, x = d14C_obs_SHZ12_ix, 
                                         y = unlist(opm_list$SHZone12)[max_ind[[2]]:6500], 
                                         z = d14C_obs[[2]]))

df_TT_SHZ12 <- data.frame(TurnoverTime = unlist(d14C_obs_SHZ12_ix),
                         Delta14C = unlist(d14C_obs[[2]]),
                         Zones = feat[2])

d14C_obs_SHZ3_ix <- map(seq_along(seq_along(d14C_obs[[3]])), 
                        ~d14C_match_fun(i = .x, x = d14C_obs_SHZ3_ix, 
                                        y = unlist(opm_list$SHZone3)[max_ind[[3]]:6500], 
                                        z = d14C_obs[[3]]))

df_TT_SHZ3 <- data.frame(TurnoverTime = unlist(d14C_obs_SHZ3_ix),
                         Delta14C = unlist(d14C_obs[[3]]),
                         Zones = feat[3])


df_TT <- rbind(df_TT_NHZ3, df_TT_SHZ12, df_TT_SHZ3)

df_14C_TT <- df_14c_zones %>% 
  left_join(df_TT, by = c("Delta14C", "Zones"))

df_14C_TT %>% 
  ggplot(aes(x = Delta14C, y = TurnoverTime)) +
  geom_point()

view(df_14C_TT)

write.csv(df_14C_TT, "./Data/df_14C_TT.csv", row.names = FALSE)


