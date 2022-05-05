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
library(bestNormalize)

#This script contains all the code to reproduce all data analysis and figures 
#in the manuscript von Fromm et al. (2022)

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
theme_own <- theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        plot.background = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = NA))

#Create facet names for depth
depth_names <- c("Topsoil" = "Topsoil (0-20 cm)",
                 "Subsoil" = "Subsoil (20-50 cm)")

#Create color variable
color_climate <- c("#ffffb2", "#c2a5cf", "#762a83", "#a6dba0", "#1b7837")

###Figure 1 and A3

##Principial component analysis
#Code inspired by:
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
#https://www.youtube.com/watch?v=Uhw-1NilmAk&list=PLnZgp6epRBbTsZEFXi_p6W48HhNyqwxIu&index=11

AfSIS_TOP <- AfSIS_14C %>%
  filter(Depth == "Topsoil")

AfSIS_BOT <- AfSIS_14C %>%
  filter(Depth == "Subsoil")

##PCA analysis
#Function to create PCA
AfSIS_PCA_fun <- function(dataset){
  dataset %>% 
    dplyr::select(CORG, GPP, Clay_2_1, Clay_1_1, Clay_8um, Mox, Quartz,
                  TurnoverTime, KG_p_group, Cultivation, Erosion) %>% 
    rename("C age" = TurnoverTime,
           "2:1 clays" = Clay_2_1,
           "C content" = CORG,
           "1:1 clays" = Clay_1_1,
           "Clay fraction" = Clay_8um,
           PCM = Mox) %>% 
    mutate_if(is.numeric, ~predict(bestNormalize::bestNormalize(.x))) %>% 
    PCA(graph = FALSE, scale.unit = FALSE, quali.sup = c(9:11), quanti.sup = 8)
}

#run function for both depths
set.seed(42)
AfSIS_PCA_TOP <- AfSIS_PCA_fun(dataset = AfSIS_TOP)
AfSIS_PCA_BOT <- AfSIS_PCA_fun(dataset = AfSIS_BOT)

#Function to extract and plot data
PCA_biplot_fun <- function(PCA_data, dataset){
  #Number of Dimensions to use
  ind <- data.frame(get_pca_ind(PCA_data)$coord)
  quanti.sup <- data.frame(PCA_data$quanti.sup$coord, row.names = "C age")
  var <- facto_summarize(PCA_data, element = "var", 
                         result = c("coord", "contrib", "cos2"))
  #factor for scaling variables to space of individuals (for Dim.1 and Dim.2)
  r <- min((max(ind[, "Dim.1"]) - min(ind[, "Dim.1"])/(max(var[, "Dim.1"]) - 
                                                         min(var[, "Dim.1"]))), 
           (max(ind[, "Dim.2"]) - min(ind[, "Dim.2"])/(max(var[, "Dim.2"]) - 
                                                         min(var[, "Dim.2"]))))
  ggplot() +
    geom_point(data = ind, aes(x = Dim.1, y = Dim.2, fill = dataset$KG_p_group),
               shape = 21, size = 7) +
    # scale factor based on fviz_pca_biplot function
    geom_segment(data = var, aes(x = 0, xend = Dim.1*r*0.7, y = 0, yend = Dim.2*r*0.7),
                 size = 1, arrow = arrow(length = unit(0.02, "npc"))) +
    geom_label(data = var, aes(x = Dim.1*r*0.75, y = Dim.2*r*0.75, label = name),
               size = 5.5, fontface = "bold", fill = "white", alpha = 0.7,
               label.size = NA, label.padding = unit(0.1, "lines")) +
    geom_segment(data = quanti.sup, aes(x = 0, xend = Dim.1*r*0.7, y = 0, yend = Dim.2*r*0.7),
                 size = 1, arrow = arrow(length = unit(0.02, "npc")),
                 color = "blue") +
    geom_label(data = quanti.sup, aes(x = Dim.1*r*0.75, y = Dim.2*r*0.75, label = "C age"),
               size = 5.5, fontface = "bold", fill = "white", alpha = 0.7,
               label.size = NA, label.padding = unit(0.1, "lines"), color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 19) +
    theme(axis.text = element_text(color = "black"),
          panel.grid =  element_blank()) 
}

#Topsoil results
AfSIS_PCA_TOP$eig
AfSIS_PCA_TOP$quali.sup$eta2

PCA_biplot_TOP <- PCA_biplot_fun(PCA_data = AfSIS_PCA_TOP, dataset = AfSIS_TOP) +
  scale_x_continuous("Dimension 1 (52.7%): Control on C content", limits = c(-4.6,4.5),
                     breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (22.5%): Control on C age") +
  # coord_fixed() +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold")) +
  scale_fill_manual("Climate zones", 
                    values = color_climate)

##Scatterplot and violin plots: Mean C age ~ SOC colored by climate zones


##Tables
