## Sophie von Fromm ##
## 2023-July-12
## Data analysis - random forest

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(mlr3spatial)
library(mlr3spatiotempcv) 
library(iml)
library(sf)
library(scales)

#Load dataset
AfSIS_14C <- read_csv("./Data/AfSIS_data_all.csv") %>% 
  dplyr::rename(Cultivation = PlotCultMgd) %>% 
  mutate(ClimateGroups = case_when(
    grepl("seasonal", KG_p_group) ~ "seasonal",
    grepl("humid", KG_p_group) ~ "humid",
    TRUE ~ "arid"
  )) %>% 
  mutate(Clay_Minerals = case_when(
    Clay_2_1 > 0 ~ "2:1 clay minerals",
    Clay_1_1 > 0 &
      Clay_2_1 == 0 ~ "1:1 clay minerals only"
  )) %>% 
  #create mineral groups for faceting
  mutate(Clay_Minerals_Content = case_when(
    Clay_Minerals == "2:1 clay minerals" ~ Clay_2_1,
    Clay_Minerals == "1:1 clay minerals only" ~ Clay_1_1
  )) %>% 
  mutate(GPP = GPP*365/1000)

AfSIS_14C$Plot <- as.character(AfSIS_14C$Plot)

AfSIS_14C$Cluster <- as.character(AfSIS_14C$Cluster)

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

climate_names <- c("Arid" = "Arid\n",
                   "Temperate (seasonal)" = "Temperate\n(seasonal)",
                   "Temperate (humid)" = "Temperate\n(humid)",
                   "Tropical (seasonal)" = "Tropical\n(seasonal)",
                   "Tropical (humid)" = "Tropical\n(humid)")

# Create empty climate zone to display legend better
AfSIS_14C$KG_p_group <- factor(AfSIS_14C$KG_p_group,
                               levels = c("Arid", "",
                                          "Temperate (seasonal)",
                                          "Tropical (seasonal)",
                                          "Temperate (humid)",
                                          "Tropical (humid)"))
##Create color variables
#Climate
color_climate <- c("#ebebcb", "#c2a5cf", "#762a83", "#a6dba0", "#1b7837")

#Create color scheme for empty climate zone
color_climate_2 <- c("#ebebcb", "white", "#c2a5cf", "#762a83", "#a6dba0", 
                     "#1b7837")

AfSIS_TOP <- AfSIS_14C %>%
  filter(Depth == "Topsoil")

AfSIS_BOT <- AfSIS_14C %>%
  filter(Depth == "Subsoil")

### Random forest analysis ###
# Filter data (only keep variables of interest)
AfSIS_rf <- AfSIS_14C %>% 
  unite("id", Country:Cluster) %>% 
  dplyr::select(id, CORG, GPP, Clay_2_1, Clay_1_1, Clay_8um, Mox, 
                Pedogenic_Oxides, Quartz, Feldspars, TurnoverTime)

task_rf <- as_task_regr(x = AfSIS_rf,
                        target = "TurnoverTime")

lrn_rf <- lrn("regr.ranger", importance = "permutation",
                  num.trees = 1000)

# Add id as group for CV (same id kept together)
task_rf$set_col_roles("id", roles = "group")
print(task_rf)

# cross-validation
set.seed(42)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(task_rf)

# autoplot(resampling_top, task_rf_top, fold_id = c(1:10)) 

## Train model & check performance
rf_all <- mlr3::resample(task = task_rf, learner = lrn_rf, 
                         resampling = resampling, store_models = TRUE)

rf_all$aggregate(measures = msrs(c("regr.rsq", "regr.mae")))

rf_pred <- rf_all$prediction(predict_sets = "test")
rf_pred_df <- data.frame(truth = rf_pred$truth,
                         response = rf_pred$response)

rf_pred_df %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(size = 1) +
  geom_rug(length = unit(0.25, "cm")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed turnover time [yr]", 
                     limits = c(-100,6200), expand = c(0,0)) +
  scale_x_continuous("Predicted turnover time [yr]", 
                     limits = c(-100,2000), expand = c(0,0))

# vi_all <- lapply(rf_all$learners, function(x) x$model$variable.importance)
# 
# vi_all %>%   
#   plyr::ldply() %>% 
#   dplyr::summarise(across(everything(), list(median = median, mad = mad), 
#                           .names = "{.col}.{.fn}")) %>% 
#   pivot_longer(everything(), names_to = "names", values_to = "values") %>% 
#   separate_wider_delim(names, ".", names = c("predictor", "measure")) %>% 
#   pivot_wider(names_from = measure, values_from = values) %>% 
#   arrange(median)

### Partial dependence plots
AfSIS_rf_pdp <- AfSIS_14C %>% 
  unite("id", Country:Cluster) %>% 
  arrange(id) %>% 
  dplyr::select(CORG, GPP, Clay_2_1, Clay_1_1, Clay_8um, Mox, 
                Pedogenic_Oxides, Quartz, Feldspars, TurnoverTime,
                id)

task_rf_pdp <- as_task_regr(x = AfSIS_rf_pdp %>% 
                              dplyr::select(-id),
                                target = "TurnoverTime")

lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                   num.trees = 1000) 

set.seed(42)
lrn_rf_pdp$train(task_rf_pdp)

model_rf <- Predictor$new(lrn_rf_pdp, data = AfSIS_rf_pdp %>% 
                            dplyr::select(-id))

#PDP
effect_rf_pdp <- FeatureEffects$new(model_rf, method = "pdp",
                                    features = c("CORG", "GPP", "Clay_2_1",
                                                 "Clay_1_1", "Clay_8um",
                                                 "Mox", "Pedogenic_Oxides",
                                                 "Quartz", "Feldspars"))

pdp_all_df <- rbindlist(effect_rf_pdp$results) %>% 
  tibble() %>% 
  dplyr::select(-.type) %>% 
  dplyr::rename(feature_value = .borders,
                predicted_value = .value,
                feature_name = .feature)

pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Clay_1_1"),
                                   "1:1 clay minerals [wt-%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Clay_2_1"),
                                   "2:1 clay minerals [wt-%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Clay_8um"),
                                   "Clay fraction [%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "CORG"),
                                   "SOC [wt-%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Feldspars"),
                                   "Feldspars [wt-%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "GPP"),
                                   "GPP [kgC/m²yr]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "GPP"),
                                   "GPP [kgC/m²yr]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Mox"),
                                   "Poorly-crystalline minerals [wt-%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Pedogenic_Oxides"),
                                   "Pedogenic oxides [wt-%]")
pdp_all_df$feature_name <- replace(pdp_all_df$feature_name,
                                   which(pdp_all_df$feature_name == "Quartz"),
                                   "Quartz [wt-%]")
  
pdp_all_df %>% 
  ggplot(aes(x = feature_value, y = predicted_value)) +
  geom_hline(yintercept = mean(AfSIS_14C$TurnoverTime), 
             color = "darkgrey", linetype = "dashed") +
  geom_line(linewidth = 1) +
  facet_wrap(~feature_name, scales = "free_x") +
  scale_y_continuous("Predicted mean SOC age [yr]", trans = reverse_trans(),
                     limits = c(1700,0), expand = c(0,0)) +
  scale_x_continuous("Range of explanatory variable", expand = c(0,0)) +
  theme_bw(base_size = 16) +
  theme_own +
  theme(axis.title = element_text(face = "plain"))

ggsave("./Figures/AfSIS_14C_FigureA9.jpeg", width = 12, height = 8)

## ICE
effect_rf_ice <- FeatureEffects$new(model_rf, method = "ice",
                                    features = c("CORG", "GPP", "Clay_2_1",
                                                 "Clay_1_1", "Clay_8um",
                                                 "Mox", "Pedogenic_Oxides",
                                                 "Quartz", "Feldspars"))

ice_all_df <- rbindlist(effect_rf_ice$results) %>% 
  tibble() %>% 
  dplyr::select(-.type) 

feature_name <- effect_rf_ice$features
names(feature_name) <- feature_name

# Function to merge model with raw data
ice_plot_df_fun <- function(x){
  AfSIS_14C %>% 
    unite("id", Country:Cluster) %>% 
    arrange(id) %>% 
    dplyr::select(id, KG_p_group) %>% 
    rownames_to_column(var = ".id") %>% 
    dplyr::mutate(.id = as.integer(.id)) %>% 
    left_join(ice_all_df[ice_all_df$.feature == x,], multiple = "all") %>% 
    dplyr::group_by(KG_p_group, .borders) %>% 
    dplyr::mutate(C_age = median(.value)) %>% 
    ungroup() %>% 
    arrange(.borders) %>% 
    dplyr::select(KG_p_group, .borders, C_age) %>% 
    distinct(.keep_all = TRUE) 
}

# Function to plot ice data
ice_plot_fun <- function(dataset){
  dataset %>% 
    ggplot(aes(x = .borders, y = C_age, color = KG_p_group)) +
    geom_path(linewidth = 2) +
    theme_bw(base_size = 19) +
    theme_own +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none",
          plot.margin = margin(l = 5, r = 0, t = 0, b = 5)) +
    scale_y_continuous("", trans = "reverse", expand = c(0,0),
                       limits = c(1800,0), breaks = seq(0,1600,400)) +
    scale_color_manual(values = color_climate)
}

line_plot_fun <- function(dataset){
  dataset %>% 
    ggplot() + 
    geom_line(aes(x = 1, y = C_age, color = KG_p_group), linewidth = 2,
              position = position_dodge(width = 1)) +
    theme_bw(base_size = 19) +
    theme_own +
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          plot.margin = margin(l = 0, r = 5, t = 0, b = 45)) +
    scale_y_continuous(trans = "reverse", expand = c(0,0),
                       limits = c(1800,0), breaks = seq(0,1600,400)) +
    scale_color_manual(values = color_climate)
}

boxplot_fun <- function(feature){
  AfSIS_14C %>% 
    ggplot(aes(x = .data[[feature]], y = KG_p_group, fill = KG_p_group)) +
    geom_boxplot() +
    theme_bw(base_size = 19) +
    theme_own +
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "top",
          plot.margin = margin(l = 66, r = 39, t = 5, b = 0),
          legend.text = element_text(size = 19))  +
    scale_fill_manual("Climate zones:", drop = FALSE,
                      values = color_climate_2) +
    guides(fill = guide_legend(ncol = 3,
                               override.aes = list(fill = color_climate_2,
                                                   color = c("black", "white",
                                                             "black", "black", 
                                                             "black", "black"))))
}

# extract legend for plotting
box_Clay_8um_plot <- boxplot_fun(feature_name["Clay_8um"]) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0))
legend_ice <- cowplot::get_legend(box_Clay_8um_plot) 

box_Quartz_plot <- boxplot_fun(feature_name["Quartz"]) +
  theme(legend.position = "right") +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_fill_manual("Climate zones:", values = color_climate) +
  guides(fill = guide_legend(ncol = 1))
legend_ice_1 <- get_legend(box_Quartz_plot)

### Arrange multiple panels together
## Minerals

# Clay content
ice_Clay_8um_p <- ice_plot_fun(ice_plot_df_fun("Clay_8um")) +
  scale_x_continuous("Clay + fine silt fraction [wt-%]", limits = c(0,100),
                     expand = c(0,0))

ice_Clay_8um_line_p <- line_plot_fun(ice_plot_df_fun("Clay_8um"))

ice_Clay_8um_range_p <- ggarrange(ice_Clay_8um_p, ice_Clay_8um_line_p, widths = c(10,1))

box_Clay_8um_plot <- boxplot_fun(feature_name["Clay_8um"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,100), expand = c(0,0))

Clay_8um_p <- ggarrange(box_Clay_8um_plot, ice_Clay_8um_range_p,  nrow = 2, 
                        heights = c(1,3))

# 2:1 clay minerals
ice_Clay_2_1_p <- ice_plot_fun(ice_plot_df_fun("Clay_2_1")) +
  scale_x_continuous("2:1 clay minerals [wt-%]", limits = c(0,65),
                     expand = c(0,0))

ice_Clay_2_1_line_p <- line_plot_fun(ice_plot_df_fun("Clay_2_1")) 

ice_Clay_2_1_range_p <- ggarrange(ice_Clay_2_1_p, ice_Clay_2_1_line_p, widths = c(10,1))

box_Clay_2_1_plot <- boxplot_fun(feature_name["Clay_2_1"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,65), expand = c(0,0))

Clay_2_1_p <- ggarrange(box_Clay_2_1_plot, ice_Clay_2_1_range_p,  nrow = 2, 
                        heights = c(1,3))

# 1:1 clay minerals
ice_Clay_1_1_p <- ice_plot_fun(ice_plot_df_fun("Clay_1_1")) +
  scale_x_continuous("1:1 clay minerals [wt-%]", limits = c(0,45),
                     expand = c(0,0)) 

ice_Clay_1_1_line_p <- line_plot_fun(ice_plot_df_fun("Clay_1_1")) 

ice_Clay_1_1_range_p <- ggarrange(ice_Clay_1_1_p, ice_Clay_1_1_line_p, widths = c(10,1))

box_Clay_1_1_plot <- boxplot_fun(feature_name["Clay_1_1"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,45), expand = c(0,0))

Clay_1_1_p <- ggarrange(box_Clay_1_1_plot, ice_Clay_1_1_range_p,  nrow = 2, 
                        heights = c(1,3))

# oxalate-extractable metals
ice_Mox_p <- ice_plot_fun(ice_plot_df_fun("Mox")) +
  scale_x_continuous("Oxalate-extractable metals [wt-%]", limits = c(0,8),
                     expand = c(0,0))

ice_Mox_line_p <- line_plot_fun(ice_plot_df_fun("Mox"))

ice_Mox_range_p <- ggarrange(ice_Mox_p, ice_Mox_line_p, widths = c(10,1))

box_Mox_plot <- boxplot_fun(feature_name["Mox"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,8), expand = c(0,0))

Mox_p <- ggarrange(box_Mox_plot, ice_Mox_range_p,  nrow = 2, 
                   heights = c(1,3))

# plot all together
ice_minerals_p <- annotate_figure(
  ggarrange(Clay_8um_p, Clay_2_1_p, Clay_1_1_p, Mox_p, 
            legend.grob = legend_ice,
            labels = c("(a)", "(b)", "(c)", "(d)")),
  left = text_grob("Predicted mean SOC age [yr]",
                   face = "bold", size = 19, rot = 90))

ggsave(ice_minerals_p, file =  "./Figures/AfSIS_14C_Figure4.jpeg", 
       width = 12, height = 8)

## C inputs
# SOC
ice_CORG_p <- ice_plot_fun(ice_plot_df_fun("CORG")) +
  scale_x_continuous("Soil organic carbon content [wt-%]", limits = c(0,8),
                     expand = c(0,0)) +
  theme_bw(base_size = 18) +
  theme_own +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "plain"),
        legend.position = "none",
        plot.margin = margin(l = 5, r = 0, t = 0, b = 5))

ice_CORG_line_p <- line_plot_fun(ice_plot_df_fun("CORG"))

ice_CORG_range_p <- ggarrange(ice_CORG_p, ice_CORG_line_p, widths = c(10,1))

box_CORG_plot <- boxplot_fun(feature_name["CORG"]) +
  scale_x_continuous(limits = c(0,8), expand = c(0,0)) +
  theme(legend.position = "none")

CORG_p <- ggarrange(box_CORG_plot, ice_CORG_range_p,  nrow = 2, 
                        heights = c(1,3))

# GPP
ice_GPP_p <- ice_plot_fun(ice_plot_df_fun("GPP")) +
  scale_x_continuous("Gross primary productivity [kgC/m²yr]", limits = c(0,2.5),
                     expand = c(0,0)) +
  theme_bw(base_size = 18) +
  theme_own +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "plain"),
        legend.position = "none",
        plot.margin = margin(l = 5, r = 0, t = 0, b = 5))

ice_GPP_line_p <- line_plot_fun(ice_plot_df_fun("GPP")) 

ice_GPP_range_p <- ggarrange(ice_GPP_p, ice_GPP_line_p, widths = c(10,1))

box_GPP_plot <- boxplot_fun(feature_name["GPP"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,2.5), expand = c(0,0))

GPP_p <- ggarrange(box_GPP_plot, ice_GPP_range_p,  nrow = 2, 
                   heights = c(1,3))

# plot both together
ice_SOC_GPP_p <- annotate_figure(
  ggarrange(CORG_p, GPP_p),
  left = text_grob("Predicted mean SOC age [yr]",
                   size = 18, rot = 90))

ggsave(ice_SOC_GPP_p, file =  "./Figures/AfSIS_14C_Figure5c.jpeg", 
       width = 12, height = 5)

## Other minerals
# Feldspars
ice_Feldspars_p <- ice_plot_fun(ice_plot_df_fun("Feldspars")) +
  scale_x_continuous("Feldspars [wt-%]", limits = c(0,55),
                     expand = c(0,0))

ice_Feldspars_line_p <- line_plot_fun(ice_plot_df_fun("Feldspars")) 

ice_Feldspars_range_p <- ggarrange(ice_Feldspars_p, ice_Feldspars_line_p, widths = c(10,1))

box_Feldspars_plot <- boxplot_fun(feature_name["Feldspars"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,55), expand = c(0,0))

Feldspars_p <- ggarrange(box_Feldspars_plot, ice_Feldspars_range_p,  nrow = 2, 
                         heights = c(1,3))

# Pedogenic oxides
ice_Pedogenic_Oxides_p <- ice_plot_fun(ice_plot_df_fun("Pedogenic_Oxides")) +
  scale_x_continuous("Pedogenic oxides [wt-%]", limits = c(0,71),
                     expand = c(0,0)) 

ice_Pedogenic_Oxides_line_p <- line_plot_fun(ice_plot_df_fun("Pedogenic_Oxides"))

ice_Pedogenic_Oxides_range_p <- ggarrange(ice_Pedogenic_Oxides_p, ice_Pedogenic_Oxides_line_p, 
                                          widths = c(10,1))

box_Pedogenic_Oxides_plot <- boxplot_fun(feature_name["Pedogenic_Oxides"]) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,71), expand = c(0,0))

Pedogenic_Oxides_p <- ggarrange(box_Pedogenic_Oxides_plot, ice_Pedogenic_Oxides_range_p,  
                                nrow = 2, heights = c(1,3))
# Quartz
ice_Quartz_p <- ice_plot_fun(ice_plot_df_fun("Quartz")) +
  scale_x_continuous("Quartz [wt-%]", limits = c(0,100),
                     expand = c(0,0)) 

ice_Quartz_line_p <- line_plot_fun(ice_plot_df_fun("Quartz")) 

ice_Quartz_range_p <- ggarrange(ice_Quartz_p, ice_Quartz_line_p, widths = c(10,1))

box_Quartz_plot <- box_Quartz_plot + 
  theme(legend.position = "none")

Quartz_p <- ggarrange(box_Quartz_plot, ice_Quartz_range_p,  nrow = 2, 
                      heights = c(1,3))

# plot all together
ice_Feld_POX_Quartz_p <- annotate_figure(
  ggarrange(Feldspars_p, Pedogenic_Oxides_p, Quartz_p, legend_ice_1,
            labels = c("a)", "b)", "c)")), 
  left = text_grob("Predicted mean SOC age [yr]",
                   face = "bold", size = 19, rot = 90))
