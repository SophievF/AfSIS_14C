## Sophie von Fromm
## 2022-Apr-08


#Function to load csv files that contain 14C data
read_csv_14C_fun <- function(path, pattern = "*.csv"){
  list.files(path, pattern, full.names = TRUE) %>% 
    map_df(~read_csv(.))
    
}




