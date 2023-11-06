## XRPD analysis for 14C samples from the AfSIS dataset##
## Benjamin Butler and Sophie von Fromm
##
## 2022-04-04 ##

library(tidyverse)

#Install the development version of powdR
devtools::install_github("benmbutler/powdR",
                         upgrade = FALSE)

library(powdR)

#Get the paths of the XRPD patterns
paths <- dir("./Data/XRPD_XY_Files", full.names = TRUE)

#Get the xy files
xrpd <- read_xy(paths)

#Only use XRPD samples that are used for 14C analysis
AfSIS_XRPD <- read_csv("./Data/AfSIS_XRPD_LongLat.csv")
AfSIS_14C <- read_csv("./Data/AfSIS_LongLat.csv") 

#Add variable names for later merging of both datasets
AfSIS_XRPD <- AfSIS_XRPD %>% 
  mutate(Site = Site.XRPD,
         Cluster = Cluster.XRPD,
         Depth = Depth.XRPD) 

#Add variable names for later merging of both datasets
#Based on 13C values we decided to exclude these 4 samples; see publication and DataPreparation.R
AfSIS_14C <- AfSIS_14C %>% 
  filter(SSN != "icr042245",
         SSN != "icr042246",
         SSN != "icr042265",
         SSN != "icr042266") %>% 
  mutate(Site.14C = Site,
         Cluster.14C = Cluster,
         Depth.14C = Depth) %>% 
  rename(SSN.14C = SSN,
         Plot.14C = Plot,
         Longitude.14C = Longitude,
         Latitude.14C = Latitude)

#Merge both datasets to see for which 14C sample we have XRPD data
AfSIS_14C_XRPD_merged <- AfSIS_14C %>% 
  left_join(AfSIS_XRPD) %>% 
  mutate(SSN.Match = SSN.14C == SSN.XRPD,
         Site.Match = Site.14C == Site.XRPD,
         Cluster.Match = Cluster.14C == Cluster.XRPD,
         Plot.Match = Plot.14C == Plot.XRPD,
         Long.Dif = abs(Longitude.14C - Longitude.XRPD),
         Lat.Dif = abs(Latitude.14C - Latitude.XRPD)) %>% 
  dplyr::select(-Site,-Cluster,-Depth)

#A few 14C samples have no matching XRPD sample (missing XRPD subsoil sample)
#Use topsoil sample instead
AfSIS_14C_XRPD_merged %>% 
  filter(is.na(SSN.XRPD))

#Fill missing value with previous row (data from topsoil layer)
AfSIS_14C_XRPD_match <- AfSIS_14C_XRPD_merged %>% 
  fill(c(SSN.XRPD:Lat.Dif)) 

#Only use unique SSN ID's for fitting
xrpd_14c <- xrpd[names(xrpd) %in% AfSIS_14C_XRPD_match$SSN.XRPD == TRUE]

#Load rockjock and afsis reference patterns from powdR
data(rockjock)
data(afsis)

rj_afsis <- merge(rockjock, afsis)

#Define the names of the minerals to include from the library
usuals <- c("Corundum", "Quartz", "K-feldspar",
            "Plagioclase", "Calcite",
            "Dolomite", "Ankerite", "Magnesite",
            "Siderite", "Halloysite", "Kaolinite",
            "Opal", "Olivine",
            "Dickite", "Smectite (Di)", "Smectite (Tri)",
            "Smectite (ML)", "Illite", "Glauconite",
            "Mica (Di)", "Mica (Tri)", "Amphibole",
            "Pyroxene", "Gypsum", "Anhydrite", "Bassanite", "Magnetite",
            "Hematite", "Goethite", "Maghemite",
            "Apatite", "Anatase", "Rutile",
            "Gibbsite", "Ilmenite", "Muscovite",
            "Cristobalite", "Palygorskite", "Vermiculite",
            "HUMIC_ACID", "FERRIHYDRITE_HUMBUG_CREEK", "FERRIHYDRITE",
            "Organic matter", "VOLCANIC_GLASS_FROM_WHITE_R_TEPHRA",
            "BACK_POS")

#subset the library
rj_afsis <- subset(rj_afsis, refs = usuals, mode = "keep")

#Define the amorphous phases
amorph <- c("HUMIC_ACID", "FERRIHYDRITE_HUMBUG_CREEK", "FERRIHYDRITE",
            "VOLCANIC_GLASS_FROM_WHITE_R_TEPHRA", "ORGANIC_MATTER",
            "ORGANIC_AFSIS")

## Parallel fitting

library(doParallel)
library(foreach)

#Detect number of cores on machine
UseCores <- detectCores()

#Register the cluster using n - 1 cores
cl <- makeCluster(UseCores-1)

registerDoParallel(cl)

#Use foreach loop and %dopar% to compute in parallel
mineral_fit <- foreach(i = 1:length(xrpd_14c)) %dopar%
  (powdR::afps(lib = rj_afsis,
               smpl = xrpd_14c[[i]],
               std = "QUARTZ_1_AFSIS",
               align = 0.2,
               lod = 0.1,
               amorphous = amorph,
               amorphous_lod = 0,
               force = "BACK_POS",
               shift = 0.05))

#name the items in the mineral_fit list
names(mineral_fit) <- names(xrpd_14c)

#stop the cluster
stopCluster(cl)

#Save the data
save(mineral_fit, 
     file = "./Data/mineral_fit.Rdata")

# Plot example spectra
plot(mineral_fit$icr033408, "Cu", interactive = TRUE)

mineral_fit$icr033408$phases
mineral_fit$icr033408$phases_grouped

#Summarize minerals and save data as csv file
mineral_grouped <- summarise_mineralogy(mineral_fit, type = "grouped", 
                                        order = TRUE)

# write.csv(mineral_grouped, row.names = FALSE,
#           "./Data/XRPD_mineral_grouped.csv")

##Merge XRPD and AfSIS data back together
#Prepare mineral data
mineral_SSN <- mineral_grouped %>% 
  dplyr::rename(SSN.XRPD = sample_id) %>% 
  replace(is.na(.), 0)

#Prepare matches between 14C and XRPD data
AfSIS_SSN_match <- AfSIS_14C_XRPD_match %>% 
  dplyr::select(SSN.14C, SSN.XRPD, Depth.14C) %>% 
  dplyr::rename(SSN = SSN.14C,
                Depth = Depth.14C) 

#Prepare 14C data
AfSIS_SSN <- read_csv("./Data/AfSIS_LongLat.csv") %>% 
  filter(SSN != "icr042245",
         SSN != "icr042246",
         SSN != "icr042265",
         SSN != "icr042266")

AfSIS_mineral_SSN <- AfSIS_SSN %>% 
  full_join(AfSIS_SSN_match, by = c("SSN", "Depth")) %>% 
  left_join(mineral_SSN, by = "SSN.XRPD")

write.csv(AfSIS_mineral_SSN, row.names = FALSE,
          "./Data/AfSIS_Mineral_fits.csv")
