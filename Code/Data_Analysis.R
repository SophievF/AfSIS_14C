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
theme_own <- theme(axis.text = element_text(color = "black"),
                   plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   strip.background = element_rect(fill = NA),
                   strip.text = element_text(face = "bold", size = 14),
                   axis.title = element_text(face = "bold"))

#Create facet names for depth and climate zones
depth_names <- c("Topsoil" = "Topsoil (0-20 cm)",
                 "Subsoil" = "Subsoil (20-50 cm)")

climate_names <- c("Arid" = "Arid\n",
                   "Temperate (seasonal)" = "Temperate\n(seasonal)",
                   "Temperate (humid)" = "Temperate\n(humid)",
                   "Tropical (seasonal)" = "Tropical\n(seasonal)",
                   "Tropical (humid)" = "Tropical\n(humid)")

##Create color variables
#Climate
color_climate <- c("#ffffb2", "#c2a5cf", "#762a83", "#a6dba0", "#1b7837")
#Erosion/Cultivation
color_Ero_Cult <- c("#80cdc1", "#dfc27d")

#Create empty climate zone to display legend better
AfSIS_14C$KG_p_group <- factor(AfSIS_14C$KG_p_group,
                               levels = c("Arid", "",
                                          "Temperate (seasonal)",
                                          "Tropical (seasonal)",
                                          "Temperate (humid)",
                                          "Tropical (humid)"))

#Create color scheme for empty climate zone
color_climate_2 <- c("#ffffb2", "white", "#c2a5cf", "#762a83", "#a6dba0", 
                     "#1b7837")

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
PCA_biplot_fun <- function(PCA_data, dataset, sub.variable){
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
    geom_point(data = ind, aes(x = Dim.1, y = Dim.2, fill = dataset),
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
    theme_own +
    theme() +
    scale_fill_manual("Climate zones", values = color_climate)
}

##Topsoil results
#Extract eigenvalues
AfSIS_PCA_TOP$eig
#Extract information about qualitative variables (climate, erosion, cultivation)
AfSIS_PCA_TOP$quali.sup
AfSIS_PCA_TOP$quali.sup$eta2
dimdesc(AfSIS_PCA_TOP, proba = 0.1)

#Plot PCA
PCA_biplot_TOP <- PCA_biplot_fun(PCA_data = AfSIS_PCA_TOP, 
                                 dataset = AfSIS_TOP$KG_p_group) +
  scale_x_continuous("Dimension 1 (52.7%): Control on SOC content", 
                     limits = c(-4.6,4.5), breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (22.5%): Control on SOC age")
  

##Subsoil results
AfSIS_PCA_BOT$eig
AfSIS_PCA_BOT$quali.sup
AfSIS_PCA_BOT$quali.sup$eta2
dimdesc(AfSIS_PCA_BOT, proba = 0.1)

PCA_biplot_BOT <- PCA_biplot_fun(PCA_data = AfSIS_PCA_BOT, 
                                 dataset = AfSIS_TOP$KG_p_group) +
  scale_x_continuous("Dimension 1 (49.5%): Control on SOC content", 
                     limits = c(-4.6,4), breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (23.1%): Control on SOC age") 

##Figure A7 and A8
#Figure A7: Cultivation
PCA_biplot_TOP_Cult <- PCA_biplot_fun(PCA_data = AfSIS_PCA_TOP, 
                                     dataset = AfSIS_TOP$Cultivation) +
  scale_x_continuous("Dimension 1 (52.7%): Control on SOC content", 
                     limits = c(-4.6,4), breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (22.5%): Control on SOC age") +
  scale_fill_manual("Cultivation", values = color_Ero_Cult)

PCA_biplot_BOT_Cult <- PCA_biplot_fun(PCA_data = AfSIS_PCA_BOT, 
                                     dataset = AfSIS_TOP$Cultivation) +
  scale_x_continuous("Dimension 1 (49.5%): Control on SOC content", 
                     limits = c(-4.6,4), breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (23.1%): Control on SOC age") +
  scale_fill_manual("Cultivation", values = color_Ero_Cult)

ggarrange(PCA_biplot_TOP_Cult, PCA_biplot_BOT_Cult, common.legend = TRUE,
          labels = c("a) Topsoil", "b) Subsoil"), label.y = 1.02)
ggsave("./Figures/AfSIS_14C_FigureA8.jpeg", width = 12, height = 8)

#Figure A8: Erosion
PCA_biplot_TOP_Ero <- PCA_biplot_fun(PCA_data = AfSIS_PCA_TOP, 
                                     dataset = AfSIS_TOP$Erosion) +
  scale_x_continuous("Dimension 1 (52.7%): Control on SOC content", 
                     limits = c(-4.6,4), breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (22.5%): Control on SOC age") +
  scale_fill_manual("Erosion", values = color_Ero_Cult)

PCA_biplot_BOT_Ero <- PCA_biplot_fun(PCA_data = AfSIS_PCA_BOT, 
                                 dataset = AfSIS_TOP$Erosion) +
  scale_x_continuous("Dimension 1 (49.5%): Control on SOC content", 
                     limits = c(-4.6,4), breaks = seq(-4,4,2)) +
  scale_y_continuous("Dimension 2 (23.1%): Control on SOC age")  +
  scale_fill_manual("Erosion", values = color_Ero_Cult)

ggarrange(PCA_biplot_TOP_Ero, PCA_biplot_BOT_Ero, common.legend = TRUE,
          labels = c("a) Topsoil", "b) Subsoil"), label.y = 1.02)
ggsave("./Figures/AfSIS_14C_FigureA8.jpeg", width = 12, height = 8)

###Scatterplot and violin plots: Mean C age ~ SOC colored by climate zones

##Scatterplot: Mean C age ~ SOC colored by climate zones
#Function for plotting
scatter_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(x = CORG, y = TurnoverTime, fill = KG_p_group)) +
    geom_point(size = 7, alpha = 0.8, shape = 21) +
    theme_bw(base_size = 19) +
    theme_own +
    theme(legend.position = "top",
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20)) +
    scale_x_continuous("SOC [wt-%]", expand = c(0,0),
                       limits = c(0,8)) +
    scale_y_continuous("Mean SOC age [yr]", trans = reverse_trans()) +
    scale_fill_manual("Climate zones:", drop = FALSE,
                      values = color_climate_2) +
    guides(fill = guide_legend(ncol = 3,
                               override.aes = list(shape = 21,
                                                   fill = color_climate_2,
                                                   color = c("black", "white",
                                                             "black", "black", 
                                                             "black", "black"))))
}

#Plot topsoil and subsoil data
scatter_TOP <- scatter_fun(dataset = AfSIS_TOP)
scatter_BOT <- scatter_fun(dataset = AfSIS_BOT)

##Violin plots
#Mean C age
vio_Cage_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(y = TurnoverTime, x = KG_p_group, fill = KG_p_group)) +
    geom_violin() +
    theme_bw(base_size = 16) +
    theme_own +
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          axis.title = element_text(face = "plain")) +
    scale_y_continuous("Mean SOC age [yr]", trans = "reverse") +
    scale_fill_manual("Climate zones",
                      values = color_climate)
}

#Plot topsoil and subsoil data
vio_Cage_TOP <- vio_Cage_fun(dataset = AfSIS_TOP)
vio_Cage_BOT <- vio_Cage_fun(dataset = AfSIS_BOT)

#Soil organic carbon
vio_SOC_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(y = CORG, x = KG_p_group, fill = KG_p_group)) +
    geom_violin() +
    theme_bw(base_size = 16) +
    theme_own +
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          axis.title = element_text(face = "plain")) +
    scale_y_continuous("SOC [wt-%]", limits = c(8,0), 
                       expand = c(0,0), trans = reverse_trans()) +
    scale_fill_manual("Climate zones",
                      values = color_climate)
}

#Plot topsoil and subsoil data
vio_SOC_TOP <- vio_SOC_fun(dataset = AfSIS_TOP)
vio_SOC_BOT <- vio_SOC_fun(dataset = AfSIS_BOT)

###Arrange plots
##Scatter and violin plots
#Topsoil
scatter_vio_TOP <- ggdraw() +
  draw_plot(scatter_TOP +
              theme(legend.position = "none")) +
  draw_plot(vio_Cage_TOP, x = 0.5, y = 0.36, width = 0.45, height = 0.28) +
  draw_plot(vio_SOC_TOP, x = 0.56, y = 0.12, width = 0.4, height = 0.25)

#Subsoil
scatter_vio_BOT <- ggdraw() +
  draw_plot(scatter_BOT +
              theme(legend.position = "none")) +
  draw_plot(vio_Cage_BOT, x = 0.5, y = 0.36, width = 0.45, height = 0.28) +
  draw_plot(vio_SOC_BOT, x = 0.56, y = 0.12, width = 0.4, height = 0.25)

##Arrange with PCA plot
#recieve legend to plot it separately
legend <- cowplot::get_legend(scatter_TOP)

#Topsoil - Figure 1
ggarrange(scatter_vio_TOP, PCA_biplot_TOP, ncol = 2, legend.grob = legend,
          labels = c("a)", "b)"), widths = c(1,1.2))
ggsave("./Figures/AfSIS_14C_Figure1.jpeg", width = 12, height = 8)

#Subsoil - Figure A3
ggarrange(scatter_vio_BOT, PCA_biplot_BOT, ncol = 2, legend.grob = legend,
          labels = c("a)", "b)"), widths = c(1,1.2))
ggsave("./Figures/AfSIS_14C_FigureA3.jpeg", width = 12, height = 8)

###Figure 2 and A4
##Scatter plot: Mean C age ~ SOC colored by Mox, clay minerals, GPP
#Mox
Mox_14C_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(y = TurnoverTime, x = CORG,
               fill = Mox)) +
    geom_point(shape = 21, size = 5) +
    facet_wrap(~KG_p_group, nrow = 1, labeller = as_labeller(climate_names)) +
    theme_bw(base_size = 20) +
    theme_own +
    theme(axis.title.x = element_blank(),
          legend.margin = margin(t = -50, r = 0, b = 0, l = 0, unit = "pt"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    scale_y_continuous("", trans = reverselog_trans(10)) +
    scale_x_continuous(trans = "log10") +
    scale_fill_viridis_c("Oxalate-extractable\nmetals [wt-%]", 
                         trans = "log10", option = "inferno", limits = c(0.01,6),
                         breaks = c(0.03,0.3,3), direction = -1) +
    guides(fill = guide_colorbar(barheight = 5, frame.colour = "black", 
                                 ticks.linewidth = 2))
}

#Plot for topsoil and subsoil samples
Mox_14C_TOP <- Mox_14C_fun(dataset = AfSIS_TOP)
Mox_14C_BOT <- Mox_14C_fun(dataset = AfSIS_BOT)

#2:1 clay minerals
Clay21_14C_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(y = TurnoverTime, x = CORG,
               fill = Clay_2_1*Clay_8um/100, color = "")) +
    geom_point(size = 5, shape = 21) +
    facet_wrap(~KG_p_group, nrow = 1) +
    theme_bw(base_size = 20) +
    theme_own +
    theme(axis.title = element_text(face = "bold"),
          strip.text = element_blank(),
          axis.title.x = element_blank(),
          legend.margin = margin(t = -15, r = 10, b = -20, l = 0, unit = "pt"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.background = element_blank()) +
    scale_y_continuous("Mean SOC age [yr]", trans = reverselog_trans(10)) +
    scale_x_continuous(trans = "log10") +
    scale_fill_viridis_c("2:1 clay minerals\nin clay fraction [%]", 
                         trans = "log10", option = "magma", direction = -1,
                         na.value = "grey80", breaks = c(0.3,3.0,30)) +
    scale_color_manual(values = "black", labels = "0") +
    guides(fill = guide_colorbar(barheight = 5, frame.colour = "black", 
                                 ticks.linewidth = 2, order = 1),
           color = guide_legend("", order = 2,
                                override.aes = list(fill = "grey80"))) 
}

#Plot for topsoil and subsoil samples
Clay21_14C_TOP <- Clay21_14C_fun(dataset = AfSIS_TOP)
Clay21_14C_BOT <- Clay21_14C_fun(dataset = AfSIS_BOT)

#GPP
GPP_14C_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(y = TurnoverTime, x = CORG,
               fill = GPP*365/1000)) +
    geom_point(shape = 21, size = 5) +
    facet_wrap(~KG_p_group, nrow = 1) +
    theme_bw(base_size = 20) +
    theme_own + 
    theme(axis.title = element_text(face = "bold"),
          strip.text = element_blank(),
          legend.margin = margin(t = 0, r = 65, b = 0, l = 0, unit = "pt"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    scale_y_continuous("", trans = reverselog_trans(10)) +
    scale_x_continuous("SOC [wt-%]", trans = "log10") +
    scale_fill_viridis_c("GPP\n[kgC/m?yr]", direction = -1, limits = c(0,2.5),
                         breaks = c(0,1,2)) +
    guides(fill = guide_colorbar(barheight = 5, frame.colour = "black", 
                                 ticks.linewidth = 2))
}

#Plot for topsoil and subsoil samples
GPP_14C_TOP <- GPP_14C_fun(dataset = AfSIS_TOP)
GPP_14C_BOT <- GPP_14C_fun(dataset = AfSIS_BOT)

ggarrange(Mox_14C_TOP, Clay21_14C_TOP, GPP_14C_TOP, nrow = 3,
          labels = c("a)", "b)", "c)"), label.x = 0.04, label.y = 1.04,
          heights = c(1.2,1,1.2))

ggsave("./Figures/AfSIS_14C_Figure2.jpeg", width = 12, height = 7)

ggarrange(Mox_14C_BOT, Clay21_14C_BOT, GPP_14C_BOT, nrow = 3,
          labels = c("a)", "b)", "c)"), label.x = 0.04, label.y = 1.04,
          heights = c(1.2,1,1.2))

ggsave("./Figures/AfSIS_14C_FigureA4.jpeg", width = 12, height = 7)

###Figure 3 and A5
##Mean C age ~ clay content colored by clay mineral type
ClayType_14C_fun <- function(dataset){
  dataset %>% 
    mutate(Clay_Minerals = case_when(
      Clay_2_1*Clay_8um/100 > 0 ~ "2:1 clay minerals",
      Clay_1_1*Clay_8um/100 > 0 &
        Clay_2_1*Clay_8um/100 == 0 ~ "1:1 clay minerals only"
    )) %>% 
    #drop samples that don't have any clay minerls
    drop_na(Clay_Minerals) %>% 
    #create mineral groups for faceting
    mutate(Clay_Minerals_Content = case_when(
      Clay_Minerals == "2:1 clay minerals" ~ Clay_2_1*Clay_8um/100,
      Clay_Minerals == "1:1 clay minerals only" ~ Clay_1_1*Clay_8um/100
    )) %>% 
    #exclude arid samples
    filter(KG_p_group != "Arid") %>% 
    ggplot(aes(y = TurnoverTime, x = Clay_8um, fill = Clay_Minerals_Content)) +
    geom_point(size = 7, shape = 21) +
    facet_wrap(~Clay_Minerals) +
    theme_bw(base_size = 21) +
    theme_own +
    theme(strip.text = element_text(size = 16),
          panel.spacing.x = unit(1, "cm"),
          legend.title = element_text(size = 19)) +
    scale_y_continuous("Mean SOC age [yr]", trans = reverselog_trans(10)) +
    scale_x_continuous("Clay + fine silt fraction [%]", expand = c(0,0),
                       limits = c(0,100)) +
    scale_fill_viridis_c("Clay mineral\ncontent in clay\nfraction [%]", 
                         option = "magma", trans = "log10", direction = -1) +
    guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                                 ticks.linewidth = 2, title.vjust = 2))
}

#Plot for topsoil and subsoil samples
ClayType_14C_fun(dataset = AfSIS_TOP)
ggsave("./AfSIS_14C_Figure3.jpeg", width = 11, height = 7)

ClayType_14C_fun(dataset = AfSIS_BOT)
ggsave("AfSIS_14C_FigureA5.jpeg", width = 11, height = 7)

###Figure A2
##Mean C age by depth

boxplot_14C <- AfSIS_14C %>%
  ggplot(aes(y = TurnoverTime, x = Depth)) +
  geom_segment(aes(x = 0.4, y = 182, xend = 1, yend = 182),
               size = 1, color = "darkgrey", linetype = "dashed") +
  geom_segment(aes(x = 0.4, y = 563, xend = 2, yend = 563),
               size = 1, color = "darkgrey", linetype = "dashed") +
  geom_boxplot(notch = TRUE) +
  theme_bw(base_size = 17) +
  theme_own +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = c("plain", "bold", "plain", "bold",
                                            "plain", "plain")),
        axis.title = element_text(face = "bold")) +
  scale_x_discrete("", labels = c("Topsoil\n(0-20 cm)", "Subsoil\n(20-50 cm)")) +
  scale_y_continuous("Mean SOC age [yr]", trans = reverselog_trans(10),
                     breaks = c(100,182,300,563,1000,3000))

boxplot_14C_diff <- AfSIS_14C %>% 
  dplyr::select(SiteID, Depth, TurnoverTime, GPP, KG_p_group) %>% 
  pivot_wider(names_from = Depth, values_from = TurnoverTime) %>% 
  mutate(TT_Dif = Topsoil-Subsoil) %>% 
  ggplot(aes(y = TT_Dif, x = KG_p_group, fill = GPP*365/1000)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  geom_jitter(size = 4, shape = 21,
              position = position_jitterdodge(0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 17) +
  theme_own +
  scale_x_discrete("", labels = climate_names) +
  scale_y_continuous("Mean SOC age [yr]: Topsoil - Subsoil", limits = c(-5200,2000),
                     expand = c(0,0), breaks = seq(-5000,2000,1000)) +
  scale_fill_viridis_c("GPP\n[kgC/m²yr]", limits = c(0,2.5), direction = -1) +
  guides(fill = guide_colorbar(barheight = 10, frame.colour = "black", 
                               ticks.linewidth = 2, title.vjust = 2))

ggarrange(boxplot_14C, boxplot_14C_diff, widths = c(0.9,2.1))

ggsave("./Figures/AfSIS_14C_FigureA2.jpeg", width = 12, height = 6)

###Figure A6
##Mean C age erosion/cultivation
#Only Topsoil data is shown which should have the strongest effect
AfSIS_TOP_EroCult <- AfSIS_TOP %>% 
  pivot_longer(cols = c(Erosion, Cultivation), values_to = "EroCult_val",
               names_to = "EroCult_name")
  
boxplot_EroCult <- AfSIS_TOP_EroCult %>%
  ggplot(aes(y = TurnoverTime, x = EroCult_name, fill = EroCult_val)) +
  facet_wrap(~Depth, scale = "free_y", labeller = as_labeller(depth_names)) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.2),
              shape = 21, size = 3) +
  theme_bw(base_size = 18) +
  theme_own + 
  scale_y_continuous("Mean SOC age [yr]", trans = "reverse") +
  scale_x_discrete("") +
  scale_fill_manual("Erosion/Cultivation", values = color_Ero_Cult)


boxplot_Cult <- AfSIS_TOP_EroCult %>%
  filter(EroCult_name == "Cultivation") %>% 
  #only use sites that have at least 3 non-cultivated and 3 cultivated plots
  filter(Site == "Analavory"|
           Site == "Didy"|
           Site == "Namasuba"|
           Site == "Nyalagari"|
           Site == "Dambidolo"|
           Site == "Pampaida") %>% 
  ggplot(aes(y = TurnoverTime, x = Site, color = EroCult_val)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.2),
              shape = 21, size = 3) +
  facet_wrap(~EroCult_name) +
  theme_bw(base_size = 18) +
  theme_own +
  scale_y_continuous("", trans = "reverse") +
  scale_x_discrete("", labels = c("Analavory,\nMadagascar",
                                  "Dambidolo,\nEthiopia",
                                  "Didy,\nMadagascar",
                                  "Namasuba,\nUganda",
                                  "Nyalagari,\nNiger",
                                  "Pampaida,\nNigeria")) +
  scale_color_manual("Erosion/Cultivation", values = color_Ero_Cult)

boxplot_Ero <- AfSIS_TOP_EroCult %>%
  filter(EroCult_name == "Erosion") %>% 
  #only use sites that have at least 3 non-eroded and 3 eroded plots
  filter(Site == "Analavory"|
           Site == "Namasuba"|
           Site == "Hoima"|
           Site == "Mucope"|
           Site == "Didy") %>% 
  ggplot(aes(y = TurnoverTime, x = Site, color = EroCult_val)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.2),
              shape = 21, size = 3) +
  facet_wrap(~EroCult_name) +
  theme_bw(base_size = 18) +
  theme_own +
  scale_y_continuous("", trans = "reverse") +
  scale_x_discrete("", labels = c("Analavory,\nMadagascar",
                                  "Didy,\nMadagascar",
                                  "Hoima,\nUganda",
                                  "Mucope,\nAngola",
                                  "Namasuba,\nUganda")) +
  scale_color_manual("Erosion/Cultivation", values = color_Ero_Cult)

ggarrange(boxplot_EroCult, ggarrange(boxplot_Cult, boxplot_Ero, nrow = 2,
                                     legend = "none"), 
          common.legend = TRUE, widths = c(1,1.6))

ggsave("./Figures/AfSIS_14C_FigureA6.jpeg", width = 12, height = 7)

##Tables
