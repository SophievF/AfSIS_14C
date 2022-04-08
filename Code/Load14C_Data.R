## Sophie von Fromm
## Merge 14C files from 14C lab
## 2022-Apr-08

library(tidyverse)

#Load csv files that contain 14C data
read_csv_fun <- function(path, pattern = "*.csv"){
  list.files(path, pattern, full.names = TRUE) %>% 
    map_df(~read_csv(.))
    
}

tbl_14c <- read_csv_fun(path = "./Data/Radiocarbon")




