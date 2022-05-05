## Sophie von Fromm ##
## 2022-May-04
## Data analysis

library(tidyverse)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(cowplot)
library(FactoMineR)
library(factoextra)

#This script contains all the code to reproduce all data analysis and figures 
#in the manuscript von Fromm et al. (2022)

#Load dataset
AfSIS_14C <- read_csv("./Data/AfSIS_data_all.csv")


AfSIS_14C$KG_p_group <- factor(AfSIS_14C$KG_p_group,
                               levels = c("Arid", 
                                          "Temperate (seasonal)",
                                          "Tropical (seasonal)",
                                          "Temperate (humid)",
                                          "Tropical (humid)"))

AfSIS_14C$Depth <- factor(AfSIS_14C$Depth, 
                          levels = c("Topsoil", "Subsoil"))

