## Gross primary productivity (GPP) data for the AfSIS data ##
## Sophie von Fromm ##
## 2022-04-27 ##

# GPP data was extracted based on long/lat from the Fluxcom network
# FLUXCOM_RS_V006	Gross Primary Productivity in	gC m-2 d-1

# Load data
AfSIS_GPP <- read_csv("./Data/GPP/GPP.Fluxcom_RS_V006.AfSIS_RefData.csv")

AfSIS_GPP_sum <- AfSIS_GPP %>% 
  #only keep values between 2001 and 2012 since samples were taken between 2009-2012
  dplyr::select(c(SiteID:`2012/12`)) %>% 
  pivot_longer(cols = c(`2001/1`:`2012/12`), names_to = "year_month",
               values_to = "GPP") %>% 
  group_by(SiteID) %>% 
  summarize(Longitude = mean(Longitude),
            Latitude = mean(Latitude),
            GPP = mean(GPP))

write.csv(AfSIS_GPP_sum, row.names = FALSE,
          "./Data/AfSIS_GPP.csv")
         