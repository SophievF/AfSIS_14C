## Sophie von Fromm
## 2022-November-02

library(tidyverse)
library(prospectr)
library(furrr)
library(FactoMineR)

#This script contains all the code to reproduce all data analysis and figures 
#in the manuscript von Fromm et al. (2023)

# all afsis samples
csv_list <- list.files(
  #this is the spectral data from Vagen et al (2020): https://data.worldagroforestry.org/dataset.xhtml?persistentId=doi:10.34725/DVN/QXCWP1
  path = "C:/Users/sfromm/Documents/GitHub/radiocarbon-afsis/data/ICRAF-AfSIS-MIR/spectra/",  
  pattern = ".csv", recursive = TRUE,
  full.names = TRUE)

spc_list <- map(csv_list, ~future_map(.x,  ~ read_csv(.x)))

spc_all <- bind_rows(spc_list) %>%
  mutate(Country = str_replace_all(Country, pattern = "Zimbambwe", replacement = "Zimbabwe")) %>%
  mutate(Country = str_replace_all(Country, pattern = "SAfrica", "South Africa"))

# afsis 14c and reference data
afsis_ref_14c <- read_csv("./Data/AfSIS_Ref_14C.csv")

afsis_ref_14c %>% 
  count(data_14c)

# Merge both datasets
spc_all_data <- spc_all %>% 
  drop_na() %>% 
  left_join(afsis_ref_14c) %>% 
  mutate(dataset = case_when(
    data_14c == TRUE ~ "AfSIS_14c",
    data_14c == FALSE ~ "AfSIS_ref",
    TRUE ~ "AfSIS_all" 
  )) %>% 
  dplyr::select(-data_14c)

spc_all_data %>% 
  dplyr::count(dataset)

# Prepare spectral data
colnames(spc_all_data) <- gsub("m", "", colnames(spc_all_data))
colnames(spc_all_data)

# extract wavenumbers
wavs <- colnames(spc_all_data)[grep("^[0-9]", colnames(spc_all_data))]

# add spectra to df
spc_all_df <- spc_all_data %>% 
  dplyr::select(SSN, Depth, Country, dataset)

spc_all_df$spc <- as.matrix(spc_all_data[,wavs])

# Pre-processing
spc_all_df$spc_pre <- spc_all_df$spc %>%
  savitzkyGolay(m = 2, p = 2, w = 17)

# plot pre-processed spectra
plot(x = as.numeric(colnames(spc_all_df$spc_pre)),
     colMeans(standardNormalVariate(spc_all_df$spc_pre)),
     col = "red", 
     type = "l", 
     xlim = c(4000, 400),
     xlab = "Wavenumber (1/cm)",
     ylab = "Pre-processed spectra")
grid(lty = 1, col = "#8080804D")

## PCA
pca_spc <- PCA(spc_all_df$spc, graph = FALSE)

pca_spc_df <- pca_spc$ind$coord %>% 
  as.data.frame()

spc_df <- cbind(pca_spc_df, dataset = spc_all_df$dataset) %>% 
  tibble()

# Function to plot density distribution
den_fun <- function(pca_data){
  pca_data %>% 
    ggplot(aes(x = Dim.1, y = Dim.2)) +
    stat_density_2d(aes(fill = after_stat(level)), color = "white", 
                    geom = "polygon", contour = TRUE, bins = 20) +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white")) +
    facet_wrap(~dataset) +
    scale_fill_viridis_c("Density", direction = -1, labels = scales::comma)
}

spc_df$dataset <- factor(spc_df$dataset,
                         levels = c("AfSIS_14c", "AfSIS_ref", "AfSIS_all"))
 
# all samples (un-processed)
den_fun(spc_df)
ggsave("./Figures/AfSIS_14C_FigureA2.jpeg", width = 10, height = 6)
