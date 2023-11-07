## Sophie von Fromm ##
## 2023-July-12
## Data analysis - linear regression

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(sf)
library(scales)
library(nlme)
library(MuMIn)

theme_own <- theme(axis.text = element_text(color = "black"),
                   panel.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   strip.background = element_rect(fill = NA),
                   strip.text = element_text(face = "bold", size = 14),
                   axis.title = element_text(face = "bold"))

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

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
  ))

AfSIS_14C$Plot <- as.factor(AfSIS_14C$Plot)

AfSIS_14C$Cluster <- as.factor(AfSIS_14C$Cluster)

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

#### Linear-mixed effect model - depth ####
AfSIS_14C_box <- AfSIS_14C %>% 
  # Only use samples with 14C and XRPD
  # filter(SSN == SSN.XRPD) %>% 
  #add small value to avoid 0's for transformation
  dplyr::mutate(Pedogenic_Oxides = Pedogenic_Oxides + 0.00001,
                Feldspars = Feldspars + 0.00001,
                Quartz = Quartz + 0.00001,
                Clay_1_1 = Clay_1_1 + 0.00001,
                Clay_2_1 = Clay_2_1 + 0.00001) %>% 
  dplyr::select(Site, Cluster, Plot, Depth, Erosion, Cultivation, KG_p_group,
                MAP, MAT, CORG, Mox, TurnoverTime, Clay_8um, Pedogenic_Oxides,
                Feldspars, Quartz, Clay_1_1, Clay_2_1, GPP) %>% 
  dplyr::mutate_if(is.numeric, ~predict(bestNormalize::boxcox(.x)))

mbx0 <- lme(TurnoverTime ~ 1,
            random = ~1|Site/Cluster/Plot, method = "ML", data = AfSIS_14C_box)

mbx1 <- update(mbx0, ~. + MAP)
mbx2 <- update(mbx1, ~. + MAT)
mbx3 <- update(mbx2, ~. + Depth)
mbx4 <- update(mbx3, ~. + Cultivation)
mbx5 <- update(mbx4, ~. + Erosion)
mbx6 <- update(mbx5, ~. + GPP)
mbx7 <- update(mbx6, ~. + CORG)
mbx8 <- update(mbx7, ~. + Clay_2_1)
mbx9 <- update(mbx8, ~. + Mox)
mbx10 <- update(mbx9, ~. + Clay_1_1)
mbx11 <- update(mbx10, ~. + Pedogenic_Oxides)
mbx12 <- update(mbx11, ~. + Clay_8um)
mbx13 <- update(mbx12, ~. + Quartz)
mbx14 <- update(mbx13, ~. + Feldspars)

#no-autocorrelation
as.data.frame(car::vif(mbx14))
as.data.frame(car::vif(mbx14)) %>% 
  filter(car::vif(mbx14) > 3)

## Summary output and diagnostic plot of full model (all parameters)
summary(mbx14)
#Fitted vs residuals
plot(mbx14, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx14, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx10, abline = c(0,1),
       main = "qqnorm Plot")

## Anova output for all models (step-wise) and full model
ava_full <- anova(mbx0,mbx1,mbx2,mbx3,mbx4,mbx5,mbx6,mbx7,mbx8,mbx9,mbx10,mbx11,
                  mbx12,mbx13,mbx14)
ava_full

FixedEffects <- c("~ 1", "... + MAP", "... + MAT", "... + Depth", "... + Cultivation", 
                  "... + Erosion", "... + GPP", "... + SOC", "... + 2:1 clays", 
                  "... + Mox", "... + 1:1 clays", "... + Pedogenic oxides",
                  "... + Clay content", "... + Quartz", "... + Feldspars")

ava_full_gt <- ava_full %>% 
  cbind(FixedEffects) %>% 
  dplyr::select(FixedEffects, df, AIC, BIC, logLik, Test, L.Ratio, 'p-value') %>% 
  gt() %>% 
  fmt_number(columns = c(3:5), decimals = 0) %>% 
  fmt_number(columns = c(7), decimals = 2) %>% 
  fmt_number(columns = c(8), decimals = 4)

gtsave(ava_full_gt, filename = "./Figures/AfSIS_14C_anova_full_model.rtf")

#MAP
MAP <- r.squaredGLMM(mbx1)[1,1]
#MAT
r.squaredGLMM(mbx2)
MAT <- r.squaredGLMM(mbx2)[1,1]-r.squaredGLMM(mbx1)[1,1]
#Depth
r.squaredGLMM(mbx3)
Depth <- r.squaredGLMM(mbx3)[1,1]-r.squaredGLMM(mbx2)[1,1]
#Cultivation
r.squaredGLMM(mbx4)
Cult <- r.squaredGLMM(mbx4)[1,1]-r.squaredGLMM(mbx3)[1,1]
#Erosion
r.squaredGLMM(mbx5)
Ero <- r.squaredGLMM(mbx5)[1,1]-r.squaredGLMM(mbx4)[1,1]
#GPP
r.squaredGLMM(mbx6)
GPP <- r.squaredGLMM(mbx6)[1,1]-r.squaredGLMM(mbx5)[1,1]
#SOC
r.squaredGLMM(mbx7)
SOC <- r.squaredGLMM(mbx7)[1,1]-r.squaredGLMM(mbx6)[1,1]
#Clay_2_1
r.squaredGLMM(mbx8)
Clay_2_1 <- r.squaredGLMM(mbx8)[1,1]-r.squaredGLMM(mbx7)[1,1]
#Mox
r.squaredGLMM(mbx9)
Mox <- r.squaredGLMM(mbx9)[1,1]-r.squaredGLMM(mbx8)[1,1]
#Clay_1_1
r.squaredGLMM(mbx10)
Clay_1_1 <- r.squaredGLMM(mbx10)[1,1]-r.squaredGLMM(mbx9)[1,1]
#Pedogenic_oxides
r.squaredGLMM(mbx11)
PedoOx <- r.squaredGLMM(mbx11)[1,1]-r.squaredGLMM(mbx10)[1,1]
#Clay_8um
r.squaredGLMM(mbx12)
Clay <- r.squaredGLMM(mbx12)[1,1]-r.squaredGLMM(mbx11)[1,1]
#Quartz
r.squaredGLMM(mbx13)
Quartz <- r.squaredGLMM(mbx13)[1,1]-r.squaredGLMM(mbx12)[1,1]
#Feldspars
r.squaredGLMM(mbx14)
Feldsp <- r.squaredGLMM(mbx14)[1,1]-r.squaredGLMM(mbx13)[1,1]

R2m.all <- tibble(MAP, MAT, Depth, Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                  Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "All samples (R² = 0.66)") %>% 
  dplyr::select(Models, everything())

### Depth
AfSIS_14C_box_d <- AfSIS_14C %>% 
  #add small value to avoid 0's for transformation
  dplyr::mutate(Pedogenic_Oxides = Pedogenic_Oxides + 0.00001,
                Feldspars = Feldspars + 0.00001,
                Quartz = Quartz + 0.00001,
                Clay_1_1 = Clay_1_1 + 0.00001,
                Clay_2_1 = Clay_2_1 + 0.00001) %>% 
  dplyr::select(Site, Cluster, Plot, Depth, Erosion, Cultivation, KG_p_group,
                MAP, MAT, CORG, Mox, TurnoverTime, Clay_8um, Pedogenic_Oxides,
                Feldspars, Quartz, Clay_1_1, Clay_2_1, GPP) %>% 
  group_by(Depth) %>% 
  dplyr::mutate_if(is.numeric, ~predict(bestNormalize::boxcox(.x)))

AfSIS_14C_box_top <- AfSIS_14C_box_d %>% 
  filter(Depth == "Topsoil")

AfSIS_14C_box_bot <- AfSIS_14C_box_d %>% 
  filter(Depth == "Subsoil")

## Topsoil
mbx0.top <- lme(TurnoverTime ~ 1,
                random = ~1|Site/Cluster, method = "ML", data = AfSIS_14C_box_top)

mbx1.top <- update(mbx0.top, ~. + MAP)
mbx2.top <- update(mbx1.top, ~. + MAT)
mbx3.top <- update(mbx2.top, ~. + Cultivation)
mbx4.top <- update(mbx3.top, ~. + Erosion)
mbx5.top <- update(mbx4.top, ~. + GPP)
mbx6.top <- update(mbx5.top, ~. + CORG)
mbx7.top <- update(mbx6.top, ~. + Clay_2_1)
mbx8.top <- update(mbx7.top, ~. + Mox)
mbx9.top <- update(mbx8.top, ~. + Clay_1_1)
mbx10.top <- update(mbx9.top, ~. + Pedogenic_Oxides)
mbx11.top <- update(mbx10.top, ~. + Clay_8um)
mbx12.top <- update(mbx11.top, ~. + Quartz)
mbx13.top <- update(mbx12.top, ~. + Feldspars)

summary(mbx13.top)
#Fitted vs residuals
plot(mbx13.top, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx13.top, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx13.top, abline = c(0,1),
       main = "qqnorm Plot")

ava_top <- anova(mbx0.top,mbx1.top,mbx2.top,mbx3.top,mbx4.top,mbx5.top,mbx6.top,
                 mbx7.top,mbx8.top,mbx9.top,mbx10.top,mbx11.top,mbx12.top,mbx13.top)
ava_top

#MAP
r.squaredGLMM(mbx1.top)
MAP <- r.squaredGLMM(mbx1.top)[1,1]
#MAT
r.squaredGLMM(mbx2.top)
MAT <- r.squaredGLMM(mbx2.top)[1,1]-r.squaredGLMM(mbx1.top)[1,1]
#Cultivation
r.squaredGLMM(mbx3.top)
Cult <- r.squaredGLMM(mbx3.top)[1,1]-r.squaredGLMM(mbx2.top)[1,1]
#Erosion
r.squaredGLMM(mbx4.top)
Ero <- r.squaredGLMM(mbx4.top)[1,1]-r.squaredGLMM(mbx3.top)[1,1]
#GPP
r.squaredGLMM(mbx5.top)
GPP <- r.squaredGLMM(mbx5.top)[1,1]-r.squaredGLMM(mbx4.top)[1,1]
#SOC
r.squaredGLMM(mbx6.top)
SOC <- r.squaredGLMM(mbx6.top)[1,1]-r.squaredGLMM(mbx5.top)[1,1]
#Clay_2_1
r.squaredGLMM(mbx7.top)
Clay_2_1 <- r.squaredGLMM(mbx7.top)[1,1]-r.squaredGLMM(mbx6.top)[1,1]
#Mox
r.squaredGLMM(mbx8.top)
Mox <- r.squaredGLMM(mbx8.top)[1,1]-r.squaredGLMM(mbx7.top)[1,1]
#Clay_1_1
r.squaredGLMM(mbx9.top)
Clay_1_1 <- r.squaredGLMM(mbx9.top)[1,1]-r.squaredGLMM(mbx8.top)[1,1]
#Pedogenic oxides
r.squaredGLMM(mbx10.top)
PedoOx <- r.squaredGLMM(mbx10.top)[1,1]-r.squaredGLMM(mbx9.top)[1,1]
#Clay content
r.squaredGLMM(mbx11.top)
Clay <- r.squaredGLMM(mbx11.top)[1,1]-r.squaredGLMM(mbx10.top)[1,1]
#Quartz
r.squaredGLMM(mbx12.top)
Quartz <- r.squaredGLMM(mbx12.top)[1,1]-r.squaredGLMM(mbx11.top)[1,1]
#Feldspars
r.squaredGLMM(mbx13.top)
Feldsp <- r.squaredGLMM(mbx13.top)[1,1]-r.squaredGLMM(mbx12.top)[1,1]

Depth <- NA

R2m.top.all <- tibble(MAP, MAT, Depth, Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                      Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Topsoil - all (R² = 0.52)") %>% 
  dplyr::select(Models, everything())

## Subsoil
mbx0.bot <- lme(TurnoverTime ~ 1,
                random = ~1|Site/Cluster, method = "ML", data = AfSIS_14C_box_bot)

mbx1.bot <- update(mbx0.bot, ~. + MAP)
mbx2.bot <- update(mbx1.bot, ~. + MAT)
mbx3.bot <- update(mbx2.bot, ~. + Cultivation)
mbx4.bot <- update(mbx3.bot, ~. + Erosion)
mbx5.bot <- update(mbx4.bot, ~. + GPP)
mbx6.bot <- update(mbx5.bot, ~. + CORG)
mbx7.bot <- update(mbx6.bot, ~. + Clay_2_1)
mbx8.bot <- update(mbx7.bot, ~. + Mox)
mbx9.bot <- update(mbx8.bot, ~. + Clay_1_1)
mbx10.bot <- update(mbx9.bot, ~. + Pedogenic_Oxides)
mbx11.bot <- update(mbx10.bot, ~. + Clay_8um)
mbx12.bot <- update(mbx11.bot, ~. + Quartz)
mbx13.bot <- update(mbx12.bot, ~. + Feldspars)

summary(mbx13.bot)
#Fitted vs residuals
plot(mbx13.bot, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx13.bot, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx13.bot, abline = c(0,1),
       main = "qqnorm Plot")

ava_bot <- anova(mbx0.bot,mbx1.bot,mbx2.bot,mbx3.bot,mbx4.bot,mbx5.bot,mbx6.bot,
                 mbx7.bot,mbx8.bot,mbx9.bot,mbx10.bot,mbx11.bot,mbx12.bot,mbx13.bot)

ava_bot

FixedEffects <- c("~ 1", "... + MAP", "... + MAT", "... + Cultivation", 
                  "... + Erosion", "... + GPP", "... + SOC", "... + 2:1 clays", 
                  "... + Mox", "... + 1:1 clays", "... + Pedogenic oxides",
                  "... + Clay content", "... + Quartz", "... + Feldspars")

ava_top_bot_gt <- rbind(ava_top, ava_bot) %>% 
  cbind(FixedEffects) %>% 
  dplyr::select(FixedEffects, df, AIC, BIC, logLik, Test, L.Ratio, 'p-value') %>% 
  gt(rowname_col = "Model") %>%
  tab_row_group(label = "Subsoil (20-50 cm)", rows = 15:28) %>% 
  tab_row_group(label = "Topsoil (0-20 cm)", rows = 1:14) %>% 
  fmt_number(columns = c(3:5), decimals = 0) %>% 
  fmt_number(columns = 7, decimals = 2) %>% 
  fmt_number(columns = 8, decimals = 4)
  
gtsave(ava_top_bot_gt, filename = "./Figures/AfSIS_14C_anova_depth_model.rtf")

#MAP
r.squaredGLMM(mbx1.bot)
MAP <- r.squaredGLMM(mbx1.bot)[1,1]
#MAT
r.squaredGLMM(mbx2.bot)
MAT <- r.squaredGLMM(mbx2.bot)[1,1]-r.squaredGLMM(mbx1.bot)[1,1]
#Cultivation
r.squaredGLMM(mbx3.bot)
Cult <- r.squaredGLMM(mbx3.bot)[1,1]-r.squaredGLMM(mbx2.bot)[1,1]
#Erosion
r.squaredGLMM(mbx4.bot)
Ero <- r.squaredGLMM(mbx4.bot)[1,1]-r.squaredGLMM(mbx3.bot)[1,1]
#GPP
r.squaredGLMM(mbx5.bot)
GPP <- r.squaredGLMM(mbx5.bot)[1,1]-r.squaredGLMM(mbx4.bot)[1,1]
#SOC
r.squaredGLMM(mbx6.bot)
SOC <- r.squaredGLMM(mbx6.bot)[1,1]-r.squaredGLMM(mbx5.bot)[1,1]
#Clay_2_1
r.squaredGLMM(mbx7.bot)
Clay_2_1 <- r.squaredGLMM(mbx7.bot)[1,1]-r.squaredGLMM(mbx6.bot)[1,1]
#Mox
r.squaredGLMM(mbx8.bot)
Mox <- r.squaredGLMM(mbx8.bot)[1,1]-r.squaredGLMM(mbx7.bot)[1,1]
#Clay_1_1
r.squaredGLMM(mbx9.bot)
Clay_1_1 <- r.squaredGLMM(mbx9.bot)[1,1]-r.squaredGLMM(mbx8.bot)[1,1]
#Pedogenic oxides
r.squaredGLMM(mbx10.bot)
PedoOx <- r.squaredGLMM(mbx10.bot)[1,1]-r.squaredGLMM(mbx9.bot)[1,1]
#Clay content
r.squaredGLMM(mbx11.bot)
Clay <- r.squaredGLMM(mbx11.bot)[1,1]-r.squaredGLMM(mbx10.bot)[1,1]
#Quartz
r.squaredGLMM(mbx12.bot)
Quartz <- r.squaredGLMM(mbx12.bot)[1,1]-r.squaredGLMM(mbx11.bot)[1,1]
#Feldspars
r.squaredGLMM(mbx13.bot)
Feldsp <- r.squaredGLMM(mbx13.bot)[1,1]-r.squaredGLMM(mbx12.bot)[1,1]

R2m.bot.all <- tibble(MAP, MAT, Depth, Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                      Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Subsoil - all (R² = 0.44)") %>% 
  dplyr::select(Models, everything())

R2m.merged <- rbind(R2m.all, R2m.top.all, R2m.bot.all)

## Plot models
R2m.merged <- rbind(R2m.all, R2m.top.all, R2m.bot.all)  %>% 
  mutate(across(MAP:Feldsp, ~ ifelse(.x < 0, 0, .x))) %>% 
  pivot_longer(!Models, names_to = "predictor", values_to = "values")

R2m.merged$Models <- factor(R2m.merged$Models, 
                            levels = c(
                              "Subsoil - all (R² = 0.44)",
                              "Topsoil - all (R² = 0.52)",
                              "All samples (R² = 0.66)"
                            ))

R2m.merged$predictor <- factor(R2m.merged$predictor,
                               levels = c("MAP", "MAT", "x",
                                          "Depth", "Cult", "Ero", 
                                          "GPP", "SOC", "y",
                                          "Clay_2_1",  "Mox", "z",
                                          "Clay_1_1", "PedoOx", "Clay", 
                                          "Quartz", "Feldsp"
                               ))

color_lmm <- c("#9ecae1", "#3182bd", "white",
               "#dadaeb", "#bcbddc", "#9e9ac8", 
               "#bae4b3", "#74c476", "white",
               "#ffffd4", "#fee391", "white",
               "#fec44f", "#fe9929", "#ec7014", 
               "#cc4c02", "#8c2d04")

R2m.merged %>% 
  ggplot(aes(y = Models, x = values*100, fill = predictor)) +
  geom_bar(stat = "identity", color = "black", 
           position = position_stack(reverse = TRUE)) +
  theme_bw(base_size = 15) +
  theme_own +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.justification = "left") +
  scale_y_discrete(labels = c("Subsoil\n(R² = 0.44)", 
                              "Topsoil\n(R² = 0.52)", 
                              "All samples\n(R² = 0.66)")) +
  scale_x_continuous("Explained variation [%]", expand = c(0,0), limits = c(0,70)) +
  scale_fill_manual(values = color_lmm, drop = FALSE,
                    labels = c("MAP","MAT", "",
                               "Depth", "Cultivation", "Erosion",
                               "GPP", "SOC", "",
                               "2:1 clays",expression(paste(M[ox])), "",
                               "1:1 clays", "Pedogenic oxides", "Clay content",
                               "Quartz", "Feldspars")) +
  guides(fill = guide_legend("Predictors", label.hjust = 0, nrow = 3,
                             override.aes = list(fill = color_lmm,
                                                 color = c("black", "black", "white",
                                                           "black", "black", "black",
                                                           "black", "black", "white",
                                                           "black", "black", "white",
                                                           "black", "black", "black",
                                                           "black", "black"))))
ggsave("./Figures/AfSIS_14C_FigureA5.jpeg", width = 9, height = 5)


### Climate model - all
AfSIS_14C_box_c <- AfSIS_14C %>% 
  #add small value to avoid 0's for transformation
  dplyr::mutate(Pedogenic_Oxides = Pedogenic_Oxides + 0.00001,
                Feldspars = Feldspars + 0.00001,
                Quartz = Quartz + 0.00001,
                Clay_1_1 = Clay_1_1 + 0.00001,
                Clay_2_1 = Clay_2_1 + 0.00001) %>% 
  dplyr::select(Site, Cluster, Plot, Depth, Erosion, Cultivation, KG_p_group,
                MAP, MAT, CORG, Mox, TurnoverTime, Clay_8um, Pedogenic_Oxides,
                Feldspars, Quartz, Clay_1_1, Clay_2_1, GPP) %>% 
  group_by(KG_p_group) %>% 
  dplyr::mutate_if(is.numeric, ~predict(bestNormalize::boxcox(.x)))

AfSIS_14C_box_arid <- AfSIS_14C_box_c %>% 
  filter(KG_p_group == "Arid")
AfSIS_14C_box_temps <- AfSIS_14C_box_c %>% 
  filter(KG_p_group == "Temperate (seasonal)")
AfSIS_14C_box_trops <- AfSIS_14C_box_c %>% 
  filter(KG_p_group == "Tropical (seasonal)")
AfSIS_14C_box_temph <- AfSIS_14C_box_c %>% 
  filter(KG_p_group == "Temperate (humid)")
AfSIS_14C_box_troph <- AfSIS_14C_box_c %>% 
  filter(KG_p_group == "Tropical (humid)")

#Arid
mbx0.arid <- lme(TurnoverTime ~ 1,
                 random = ~1|Site/Cluster/Plot, method = "ML", data = AfSIS_14C_box_arid)

mbx1.arid <- update(mbx0.arid, ~. + Cultivation)
mbx2.arid <- update(mbx1.arid, ~. + Erosion)
mbx3.arid <- update(mbx2.arid, ~. + GPP)
mbx4.arid <- update(mbx3.arid, ~. + CORG)
mbx5.arid <- update(mbx4.arid, ~. + Clay_2_1)
mbx6.arid <- update(mbx5.arid, ~. + Mox)
mbx7.arid <- update(mbx6.arid, ~. + Clay_1_1)
mbx8.arid <- update(mbx7.arid, ~. + Pedogenic_Oxides)
mbx9.arid <- update(mbx8.arid, ~. + Clay_8um)
mbx10.arid <- update(mbx9.arid, ~. + Quartz)
mbx11.arid <- update(mbx10.arid, ~. + Feldspars)

summary(mbx11.arid)
#Fitted vs residuals
plot(mbx11.arid, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx11.arid, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx11.arid, abline = c(0,1),
       main = "qqnorm Plot")

ava_arid <- anova(mbx0.arid,mbx1.arid,mbx2.arid,mbx3.arid,mbx4.arid,mbx5.arid,
                  mbx6.arid,mbx7.arid,mbx8.arid,mbx9.arid,mbx10.arid,mbx11.arid)

ava_arid

#Cultivation
r.squaredGLMM(mbx1.arid)
Cult <- r.squaredGLMM(mbx1.arid)[1,1]-r.squaredGLMM(mbx0.arid)[1,1]
#Erosion
r.squaredGLMM(mbx2.arid)
Ero <- r.squaredGLMM(mbx2.arid)[1,1]-r.squaredGLMM(mbx1.arid)[1,1]
#GPP
r.squaredGLMM(mbx3.arid)
GPP <- r.squaredGLMM(mbx3.arid)[1,1]-r.squaredGLMM(mbx2.arid)[1,1]
#SOC
r.squaredGLMM(mbx4.arid)
SOC <- r.squaredGLMM(mbx4.arid)[1,1]-r.squaredGLMM(mbx3.arid)[1,1]
#Clay_2_1
r.squaredGLMM(mbx5.arid)
Clay_2_1 <- r.squaredGLMM(mbx5.arid)[1,1]-r.squaredGLMM(mbx4.arid)[1,1]
#Mox
r.squaredGLMM(mbx6.arid)
Mox <- r.squaredGLMM(mbx6.arid)[1,1]-r.squaredGLMM(mbx5.arid)[1,1]
#Clay_1_1
r.squaredGLMM(mbx7.arid)
Clay_1_1 <- r.squaredGLMM(mbx7.arid)[1,1]-r.squaredGLMM(mbx6.arid)[1,1]
#Pedogenic_oxides
r.squaredGLMM(mbx8.arid)
PedoOx <- r.squaredGLMM(mbx8.arid)[1,1]-r.squaredGLMM(mbx7.arid)[1,1]
#Clay_8um
r.squaredGLMM(mbx9.arid)
Clay <- r.squaredGLMM(mbx9.arid)[1,1]-r.squaredGLMM(mbx8.arid)[1,1]
#Quartz
r.squaredGLMM(mbx10.arid)
Quartz <- r.squaredGLMM(mbx10.arid)[1,1]-r.squaredGLMM(mbx9.arid)[1,1]
#Feldspars
r.squaredGLMM(mbx11.arid)
Feldsp <- r.squaredGLMM(mbx11.arid)[1,1]-r.squaredGLMM(mbx10.arid)[1,1]

R2m.arid <- tibble(Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                   Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Arid (R² = 0.63)") %>% 
  dplyr::select(Models, everything())

#Temperate (seasonal)
mbx0.temps <- lme(TurnoverTime ~ 1,
                 random = ~1|Site/Cluster/Plot, method = "ML", data = AfSIS_14C_box_temps)

mbx1.temps <- update(mbx0.temps, ~. + Cultivation)
mbx2.temps <- update(mbx1.temps, ~. + Erosion)
mbx3.temps <- update(mbx2.temps, ~. + GPP)
mbx4.temps <- update(mbx3.temps, ~. + CORG)
mbx5.temps <- update(mbx4.temps, ~. + Clay_2_1)
mbx6.temps <- update(mbx5.temps, ~. + Mox)
mbx7.temps <- update(mbx6.temps, ~. + Clay_1_1)
mbx8.temps <- update(mbx7.temps, ~. + Pedogenic_Oxides)
mbx9.temps <- update(mbx8.temps, ~. + Clay_8um)
mbx10.temps <- update(mbx9.temps, ~. + Quartz)
mbx11.temps <- update(mbx10.temps, ~. + Feldspars)

summary(mbx11.temps)
#Fitted vs residuals
plot(mbx11.temps, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx11.temps, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx11.temps, abline = c(0,1),
       main = "qqnorm Plot")

ava_temps <- anova(mbx0.temps,mbx1.temps,mbx2.temps,mbx3.temps,mbx4.temps,mbx5.temps,
                   mbx6.temps,mbx7.temps,mbx8.temps,mbx9.temps,mbx10.temps,mbx11.temps)
ava_temps

#Cultivation
r.squaredGLMM(mbx1.temps)
Cult <- r.squaredGLMM(mbx1.temps)[1,1]-r.squaredGLMM(mbx0.temps)[1,1]
#Erosion
r.squaredGLMM(mbx2.temps)
Ero <- r.squaredGLMM(mbx2.temps)[1,1]-r.squaredGLMM(mbx1.temps)[1,1]
#GPP
r.squaredGLMM(mbx3.temps)
GPP <- r.squaredGLMM(mbx3.temps)[1,1]-r.squaredGLMM(mbx2.temps)[1,1]
#SOC
r.squaredGLMM(mbx4.temps)
SOC <- r.squaredGLMM(mbx4.temps)[1,1]-r.squaredGLMM(mbx3.temps)[1,1]
#Clay_2_1
r.squaredGLMM(mbx5.temps)
Clay_2_1 <- r.squaredGLMM(mbx5.temps)[1,1]-r.squaredGLMM(mbx4.temps)[1,1]
#Mox
r.squaredGLMM(mbx6.temps)
Mox <- r.squaredGLMM(mbx6.temps)[1,1]-r.squaredGLMM(mbx5.temps)[1,1]
#Clay_1_1
r.squaredGLMM(mbx7.temps)
Clay_1_1 <- r.squaredGLMM(mbx7.temps)[1,1]-r.squaredGLMM(mbx6.temps)[1,1]
#Pedogenic_oxides
r.squaredGLMM(mbx8.temps)
PedoOx <- r.squaredGLMM(mbx8.temps)[1,1]-r.squaredGLMM(mbx8.temps)[1,1]
#Clay_8um
r.squaredGLMM(mbx9.temps)
Clay <- r.squaredGLMM(mbx9.temps)[1,1]-r.squaredGLMM(mbx8.temps)[1,1]
#Quartz
r.squaredGLMM(mbx10.temps)
Quartz <- r.squaredGLMM(mbx10.temps)[1,1]-r.squaredGLMM(mbx9.temps)[1,1]
#Feldspars
r.squaredGLMM(mbx11.temps)
Feldsp <- r.squaredGLMM(mbx11.temps)[1,1]-r.squaredGLMM(mbx10.temps)[1,1]

R2m.temps <- tibble(Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                    Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Temp (seas) (R² = 0.70)") %>% 
  dplyr::select(Models, everything())

#Tropical (seasonal)
mbx0.trops <- lme(TurnoverTime ~ 1,
                  random = ~1|Site/Cluster/Plot, method = "ML", data = AfSIS_14C_box_trops)

mbx1.trops <- update(mbx0.trops, ~. + Cultivation)
mbx2.trops <- update(mbx1.trops, ~. + Erosion)
mbx3.trops <- update(mbx2.trops, ~. + GPP)
mbx4.trops <- update(mbx3.trops, ~. + CORG)
mbx5.trops <- update(mbx4.trops, ~. + Clay_2_1)
mbx6.trops <- update(mbx5.trops, ~. + Mox)
mbx7.trops <- update(mbx6.trops, ~. + Clay_1_1)
mbx8.trops <- update(mbx7.trops, ~. + Pedogenic_Oxides)
mbx9.trops <- update(mbx8.trops, ~. + Clay_8um)
mbx10.trops <- update(mbx9.trops, ~. + Quartz)
mbx11.trops <- update(mbx10.trops, ~. + Feldspars)

summary(mbx11.trops)
#Fitted vs residuals
plot(mbx11.trops, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx11.trops, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx11.trops, abline = c(0,1),
       main = "qqnorm Plot")

ava_trops <- anova(mbx0.trops,mbx1.trops,mbx2.trops,mbx3.trops,mbx4.trops,mbx5.trops,
                   mbx6.trops,mbx7.trops,mbx8.trops,mbx9.trops,mbx10.trops,mbx11.trops)
ava_trops

#Cultivation
r.squaredGLMM(mbx1.trops)
Cult <- r.squaredGLMM(mbx1.trops)[1,1]-r.squaredGLMM(mbx0.trops)[1,1]
#Erosion
r.squaredGLMM(mbx2.trops)
Ero <- r.squaredGLMM(mbx2.trops)[1,1]-r.squaredGLMM(mbx1.trops)[1,1]
#GPP
r.squaredGLMM(mbx3.trops)
GPP <- r.squaredGLMM(mbx3.trops)[1,1]-r.squaredGLMM(mbx2.trops)[1,1]
#SOC
r.squaredGLMM(mbx4.trops)
SOC <- r.squaredGLMM(mbx4.trops)[1,1]-r.squaredGLMM(mbx3.trops)[1,1]
#Clay_2_1
r.squaredGLMM(mbx5.trops)
Clay_2_1 <- r.squaredGLMM(mbx5.trops)[1,1]-r.squaredGLMM(mbx4.trops)[1,1]
#Mox
r.squaredGLMM(mbx6.trops)
Mox <- r.squaredGLMM(mbx6.trops)[1,1]-r.squaredGLMM(mbx5.trops)[1,1]
#Clay_1_1
r.squaredGLMM(mbx7.trops)
Clay_1_1 <- r.squaredGLMM(mbx7.trops)[1,1]-r.squaredGLMM(mbx6.trops)[1,1]
#Pedogenic_oxides
r.squaredGLMM(mbx8.trops)
PedoOx <- r.squaredGLMM(mbx8.trops)[1,1]-r.squaredGLMM(mbx7.trops)[1,1]
#Clay_8um
r.squaredGLMM(mbx9.trops)
Clay <- r.squaredGLMM(mbx9.trops)[1,1]-r.squaredGLMM(mbx8.trops)[1,1]
#Quartz
r.squaredGLMM(mbx10.trops)
Quartz <- r.squaredGLMM(mbx10.trops)[1,1]-r.squaredGLMM(mbx9.trops)[1,1]
#Feldspars
r.squaredGLMM(mbx11.trops)
Feldsp <- r.squaredGLMM(mbx11.trops)[1,1]-r.squaredGLMM(mbx10.trops)[1,1]


R2m.trops <- tibble(Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                    Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Trop (seas) (R² = 0.65)") %>% 
  dplyr::select(Models, everything())

#Temperate (humid)
mbx0.temph <- lme(TurnoverTime ~ 1,
                  random = ~1|Site/Cluster/Plot, method = "ML", data = AfSIS_14C_box_temph)

mbx1.temph <- update(mbx0.temph, ~. + Cultivation)
mbx2.temph <- update(mbx1.temph, ~. + Erosion)
mbx3.temph <- update(mbx2.temph, ~. + GPP)
mbx4.temph <- update(mbx3.temph, ~. + CORG)
mbx5.temph <- update(mbx4.temph, ~. + Clay_2_1)
mbx6.temph <- update(mbx5.temph, ~. + Mox)
mbx7.temph <- update(mbx6.temph, ~. + Clay_1_1)
mbx8.temph <- update(mbx7.temph, ~. + Pedogenic_Oxides)
mbx9.temph <- update(mbx8.temph, ~. + Clay_8um)
mbx10.temph <- update(mbx9.temph, ~. + Quartz)
mbx11.temph <- update(mbx10.temph, ~. + Feldspars)

summary(mbx11.temph)
#Fitted vs residuals
plot(mbx11.temph, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx11.temph, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx11.temph, abline = c(0,1),
       main = "qqnorm Plot")

ava_temph <- anova(mbx0.temph,mbx1.temph,mbx2.temph,mbx3.temph,mbx4.temph,mbx5.temph,
                   mbx6.temph,mbx7.temph,mbx8.temph,mbx9.temph,mbx10.temph,mbx11.temph)
ava_temph

#Cultivation
r.squaredGLMM(mbx1.temph)
Cult <- r.squaredGLMM(mbx1.temph)[1,1]-r.squaredGLMM(mbx2.temph)[1,1]
#Erosion
r.squaredGLMM(mbx2.temph)
Ero <- r.squaredGLMM(mbx2.temph)[1,1]-r.squaredGLMM(mbx1.temph)[1,1]
#GPP
r.squaredGLMM(mbx3.temph)
GPP <- r.squaredGLMM(mbx3.temph)[1,1]-r.squaredGLMM(mbx2.temph)[1,1]
#SOC
r.squaredGLMM(mbx4.temph)
SOC <- r.squaredGLMM(mbx4.temph)[1,1]-r.squaredGLMM(mbx3.temph)[1,1]
#Clay_2_1
r.squaredGLMM(mbx5.temph)
Clay_2_1 <- r.squaredGLMM(mbx5.temph)[1,1]-r.squaredGLMM(mbx4.temph)[1,1]
#Mox
r.squaredGLMM(mbx6.temph)
Mox <- r.squaredGLMM(mbx6.temph)[1,1]-r.squaredGLMM(mbx5.temph)[1,1]
#Clay_1_1
r.squaredGLMM(mbx7.temph)
Clay_1_1 <- r.squaredGLMM(mbx7.temph)[1,1]-r.squaredGLMM(mbx6.temph)[1,1]
#Pedogenic_oxides
r.squaredGLMM(mbx8.temph)
PedoOx <- r.squaredGLMM(mbx8.temph)[1,1]-r.squaredGLMM(mbx7.temph)[1,1]
#Clay_8um
r.squaredGLMM(mbx9.temph)
Clay <- r.squaredGLMM(mbx9.temph)[1,1]-r.squaredGLMM(mbx8.temph)[1,1]
#Quartz
r.squaredGLMM(mbx10.temph)
Quartz <- r.squaredGLMM(mbx10.temph)[1,1]-r.squaredGLMM(mbx9.temph)[1,1]
#Feldspars
r.squaredGLMM(mbx11.temph)
Feldsp <- r.squaredGLMM(mbx11.temph)[1,1]-r.squaredGLMM(mbx10.temph)[1,1]


R2m.temph <- tibble(Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                    Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Temp (humid) (R² = 0.70)") %>% 
  dplyr::select(Models, everything())

#Tropical (humid)
mbx0.troph <- lme(TurnoverTime ~ 1,
                  random = ~1|Site/Cluster/Plot, method = "ML", data = AfSIS_14C_box_troph)

mbx1.troph <- update(mbx0.troph, ~. + Cultivation)
mbx2.troph <- update(mbx1.troph, ~. + Erosion)
mbx3.troph <- update(mbx2.troph, ~. + GPP)
mbx4.troph <- update(mbx3.troph, ~. + CORG)
mbx5.troph <- update(mbx4.troph, ~. + Clay_2_1)
mbx6.troph <- update(mbx5.troph, ~. + Mox)
mbx7.troph <- update(mbx6.troph, ~. + Clay_1_1)
mbx8.troph <- update(mbx7.troph, ~. + Pedogenic_Oxides)
mbx9.troph <- update(mbx8.troph, ~. + Clay_8um)
mbx10.troph <- update(mbx9.troph, ~. + Quartz)
mbx11.troph <- update(mbx10.troph, ~. + Feldspars)

summary(mbx11.troph)
#Fitted vs residuals
plot(mbx11.troph, main = "Residuals vs Fitted")
#scale-location plot
plot(mbx11.troph, sqrt(abs(resid(.))) ~ fitted(.))
#Q-Q-Plot
qqnorm(mbx11.troph, abline = c(0,1),
       main = "qqnorm Plot")

ava_troph <- anova(mbx0.troph,mbx1.troph,mbx2.troph,mbx3.troph,mbx4.troph,mbx5.troph,
                   mbx6.troph,mbx7.troph,mbx8.troph,mbx9.troph,mbx10.troph,mbx11.troph)

FixedEffects <- c("~ 1", "... + Cultivation", "... + Erosion", "... + GPP", 
                  "... + SOC", "... + 2:1 clays",  "... + Mox", "... + 1:1 clays", 
                  "... + Pedogenic oxides", "... + Clay content", "... + Quartz", 
                  "... + Feldspars")

ava_climate_gt <- rbind(ava_arid, ava_temps, ava_trops, ava_temph, ava_troph) %>% 
  cbind(FixedEffects) %>% 
  dplyr::select(FixedEffects, df, AIC, BIC, logLik, Test, L.Ratio, 'p-value') %>% 
  gt(rowname_col = "Model") %>%
  tab_row_group(label = "Tropical (humid)", rows = 49:60) %>%
  tab_row_group(label = "Temperate (humid)", rows = 37:48) %>%
  tab_row_group(label = "Tropical (seasonal)", rows = 25:36) %>%
  tab_row_group(label = "Temperate (seasonal)", rows = 13:24) %>% 
  tab_row_group(label = "Arid", rows = 1:12) %>% 
  fmt_number(columns = c(3:5), decimals = 0) %>% 
  fmt_number(columns = 7, decimals = 2) %>% 
  fmt_number(columns = 8, decimals = 4)

gtsave(ava_climate_gt, filename = "./Figures/AfSIS_14C_anova_climate_model.rtf")

#Cultivation
r.squaredGLMM(mbx1.troph)
Cult <- r.squaredGLMM(mbx1.troph)[1,1]-r.squaredGLMM(mbx0.troph)[1,1]
#Erosion
r.squaredGLMM(mbx2.troph)
Ero <- r.squaredGLMM(mbx2.troph)[1,1]-r.squaredGLMM(mbx1.troph)[1,1]
#GPP
r.squaredGLMM(mbx3.troph)
GPP <- r.squaredGLMM(mbx3.troph)[1,1]-r.squaredGLMM(mbx2.troph)[1,1]
#SOC
r.squaredGLMM(mbx4.troph)
SOC <- r.squaredGLMM(mbx4.troph)[1,1]-r.squaredGLMM(mbx3.troph)[1,1]
#Clay_2_1
r.squaredGLMM(mbx5.troph)
Clay_2_1 <- r.squaredGLMM(mbx5.troph)[1,1]-r.squaredGLMM(mbx4.troph)[1,1]
#Mox
r.squaredGLMM(mbx6.troph)
Mox <- r.squaredGLMM(mbx6.troph)[1,1]-r.squaredGLMM(mbx5.troph)[1,1]
#Clay_1_1
r.squaredGLMM(mbx7.troph)
Clay_1_1 <- r.squaredGLMM(mbx7.troph)[1,1]-r.squaredGLMM(mbx6.troph)[1,1]
#Pedogenic_oxides
r.squaredGLMM(mbx8.troph)
PedoOx <- r.squaredGLMM(mbx8.troph)[1,1]-r.squaredGLMM(mbx7.troph)[1,1]
#Clay_8um
r.squaredGLMM(mbx9.troph)
Clay <- r.squaredGLMM(mbx9.troph)[1,1]-r.squaredGLMM(mbx8.troph)[1,1]
#Quartz
r.squaredGLMM(mbx10.troph)
Quartz <- r.squaredGLMM(mbx10.troph)[1,1]-r.squaredGLMM(mbx9.troph)[1,1]
#Feldspars
r.squaredGLMM(mbx11.troph)
Feldsp <- r.squaredGLMM(mbx11.troph)[1,1]-r.squaredGLMM(mbx10.troph)[1,1]

R2m.troph <- tibble(Cult, Ero, GPP, SOC, Clay_2_1, Mox,
                    Clay_1_1, PedoOx, Clay, Quartz, Feldsp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Trop (humid) (R² = 0.77)") %>% 
  dplyr::select(Models, everything())

## Plot all models together (Figure 3 and A8)
R2m.merged <- rbind(R2m.arid, R2m.temps, R2m.trops, R2m.temph, R2m.troph)  %>% 
  mutate(across(Cult:Feldsp, ~ ifelse(.x < 0, 0, .x))) %>% 
  pivot_longer(!Models, names_to = "predictor", values_to = "values")

R2m.merged$Models <- factor(R2m.merged$Models, 
                            levels = c(
                              "Trop (humid) (R² = 0.77)",
                              "Temp (humid) (R² = 0.70)",
                              "Trop (seas) (R² = 0.65)",
                              "Temp (seas) (R² = 0.70)",
                              "Arid (R² = 0.63)"
                            ))

R2m.merged$predictor <- factor(R2m.merged$predictor,
                               levels = c("Cult", "Ero", "GPP", "SOC", "Clay_2_1",
                                          "Mox", "Clay_1_1", "PedoOx", "Clay", 
                                          "Quartz", "Feldsp"
                               ))

color_lmm <- c("#bcbddc", "#9e9ac8", "#bae4b3", "#74c476",
               "#ffffd4", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02",
               "#8c2d04")

p_LMM_climate <- R2m.merged %>% 
  ggplot(aes(y = Models, x = values*100, fill = predictor)) +
  geom_bar(stat = "identity", color = "black",
           position = position_stack(reverse = TRUE)) +
  theme_bw(base_size = 16) +
  theme_own +
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        legend.background = element_blank(),
        axis.text.y = element_text(face = "bold"),
        legend.margin = margin(0.7,0,0,0, "cm")) +
  scale_y_discrete(labels = c("Trop (humid)\nR² = 0.77",
                              "Temp (humid)\nR² = 0.70",
                              "Trop (seas)\nR² = 0.65",
                              "Temp (seas)\nR² = 0.70",
                              "Arid\nR² = 0.63")) +
  scale_x_continuous("Explained variation [%]", expand = c(0,0), limits = c(0,80)) +
  scale_fill_manual(values = color_lmm,
                    labels = c("Cultivation", "Erosion",
                               "GPP", "SOC",
                               "2:1 clays", expression(paste(M[ox])),
                               "1:1 clays", "Pedogenic oxides",
                               "Clay content", "Quartz",
                               "Feldspars")) +
  guides(fill = guide_legend("Predictors", direction = "horizontal", nrow = 2,
                             label.hjust = 0))
# ggsave("./Figures/AfSIS_14C_LMM_Climate_bar.jpeg", width = 11, height = 6)

p_GPP_climate <- AfSIS_14C %>% 
  filter(Depth == "Topsoil") %>% 
  ggplot(aes(x = CORG, y = TurnoverTime, fill = GPP*365/1000)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme_own +
  theme(legend.position = "top",
        legend.box = "vertical", 
        legend.margin = margin(1,0,0,0, "cm"),
        strip.text = element_text(size = 10, face = ("plain"))) +
  facet_wrap(~KG_p_group, ncol = 1) +
  scale_y_continuous("Mean soil organic carbon age [yr]", 
                     trans = reverselog_trans(), breaks = c(100,300,1000,3000)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_fill_viridis_c("GPP [kgC/m²yr]", direction = -1, limits = c(0,2.6),
                       breaks = c(0,1,2)) +
  guides(fill = guide_colorbar(barheight = 1, frame.colour = "black", 
                               ticks.linewidth = 2, order = 1,
                               title.vjust = 1))

p_Mox_climate <- AfSIS_14C %>% 
  filter(Depth == "Topsoil") %>%
  ggplot(aes(x = CORG, y = TurnoverTime, fill = Mox)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme_own +
  theme(legend.position = "top",
        legend.box = "vertical", 
        legend.margin = margin(1,0,0,0, "cm"),
        strip.text = element_text(size = 10, face = ("plain"))) +
  facet_wrap(~KG_p_group, ncol = 1) +
  scale_y_continuous("", trans = reverselog_trans(), breaks = c(100,300,1000,3000)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_fill_viridis_c(expression(paste(M[ox], " [wt-%]")), 
                       trans = "log10", option = "inferno", limits = c(0.01,10),
                       breaks = c(0.01,0.1,1,10), direction = -1) +
  guides(fill = guide_colorbar(barheight = 1, frame.colour = "black", 
                               ticks.linewidth = 2, order = 1,
                               title.vjust = 1))

p_GPP_Mox_climate <- ggarrange(p_GPP_climate, p_Mox_climate)

ggarrange(p_LMM_climate, p_GPP_Mox_climate, widths = c(1.6,1), 
          labels = c("a) Linear-mixed effects models - all samples",
                     "b) Raw data - topsoil samples only"), label.x = c(-0.15,-0.2))

ggsave("./Figures/AfSIS_14C_Figure3.jpeg", width = 17.6, height = 11)

## Subsoil data only (Figure A8)
p_GPP_climate_bot <- AfSIS_14C %>% 
  filter(Depth == "Subsoil") %>% 
  ggplot(aes(x = CORG, y = TurnoverTime, fill = GPP*365/1000)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme_own +
  theme(legend.position = "top",
        legend.box = "vertical", 
        legend.margin = margin(1,0,0,0, "cm"),
        strip.text = element_text(size = 10, face = ("plain"))) +
  facet_wrap(~KG_p_group, ncol = 1) +
  scale_y_continuous("Mean soil organic carbon age [yr]", 
                     trans = reverselog_trans(), breaks = c(100,300,1000,3000)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_fill_viridis_c("GPP [kgC/m²yr]", direction = -1, limits = c(0,2.6),
                       breaks = c(0,1,2)) +
  guides(fill = guide_colorbar(barheight = 1, frame.colour = "black", 
                               ticks.linewidth = 2, order = 1,
                               title.vjust = 1))

p_Mox_climate_bot <- AfSIS_14C %>% 
  filter(Depth == "Subsoil") %>%
  ggplot(aes(x = CORG, y = TurnoverTime, fill = Mox)) +
  geom_point(size = 3, shape = 21) +
  theme_bw(base_size = 16) +
  theme_own +
  theme(legend.position = "top",
        legend.box = "vertical", 
        legend.margin = margin(1,0,0,0, "cm"),
        strip.text = element_text(size = 10, face = ("plain"))) +
  facet_wrap(~KG_p_group, ncol = 1) +
  scale_y_continuous("", trans = reverselog_trans(), breaks = c(100,300,1000,3000)) +
  scale_x_continuous("SOC [wt-%]", trans = "log10") +
  scale_fill_viridis_c(expression(paste(M[ox], " [wt-%]")), 
                       trans = "log10", option = "inferno", limits = c(0.01,10),
                       breaks = c(0.01,0.1,1,10), direction = -1) +
  guides(fill = guide_colorbar(barheight = 1, frame.colour = "black", 
                               ticks.linewidth = 2, order = 1,
                               title.vjust = 1))
ggarrange(p_GPP_climate_bot, p_Mox_climate_bot)
ggsave("./Figures/AfSIS_14C_FigureA8.jpeg", width = 8, height = 10)
