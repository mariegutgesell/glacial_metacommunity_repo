##Create species matrix arrays for fish
##M.Gutgesell
##Dec 12, 2024

library(tidyverse)
library(readxl)

fish_df <- read_excel("data/SEAK_Fish_Biomass_Supplement.xlsx")



##Add month_date to standardize to sampling month 
fish_df$Date <- as.Date(fish_df$Date, "%Y-%m-%d")
fish_df <- fish_df %>% mutate(Month_Year = format(Date, "%Y-%m")) %>%
  filter(Stream != "Transitional") %>%
  filter(Month_Year < "2019-01") %>%
  mutate(species2 = str_extract(Species, "^(.*?)(?=Juvenile)")) %>% 
  mutate(species3 = str_extract(Species, "^(.*?)(?=Fry)")) %>% 
  mutate(species2 = coalesce(species2, species3)) %>%
  mutate(species2 = coalesce(species2, Species)) %>%
  separate(`95% CI`, into = c("lower_95_CI", "upper_95_CI"), sep = "-")


fish_matrix <- xtabs(Biomass ~ species2 + Month_Year + Stream, data = fish_df)

saveRDS(fish_matrix, "data/intermediate_data/fish_sp_matrix.rds")


##Create a list of arrays with all unique possible site combinations (not including mixed site)-- Fish --------
sites <- unique(fish_df$Stream)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")
# Initialize an empty list to store the xtabs arrays
fish_xtabs_list <- list()

#Loop through each combination and create an xtabs array
for (i in 1:nrow(unique_combinations)) {
  combo <- unique_combinations[i, ]
  selected_sites <- c(combo$site1, combo$site2, combo$site3)
  
  site_counts <- table(selected_sites)
  
  replicated_data <- data.frame()
  
  #For each site in the combination, replicate its data based on how many times it occurs
  for (site in names(site_counts)) {
    count <- site_counts[[site]]
    site_data <- fish_df[fish_df$Stream == site, ]
    
    # Replicate the data according to how many times the site appears in the combination
    for (replicate_number in 1:count) {
      
      replicated_site_data <- site_data
      replicated_site_data$Replicate_Type <- replicate_number
      
      replicated_data <- rbind(replicated_data, replicated_site_data) 
    }
    
  }
  replicated_data <- replicated_data %>%
    mutate(streamID_2 = paste(Stream, Replicate_Type, sep = "_"))
  #Create an xtabs array for the replicated dataframe
  xtabs_array <- xtabs(Biomass ~ Species + Month_Year + streamID_2, data = replicated_data)
  
  fish_xtabs_list[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
}

saveRDS(fish_xtabs_list, "data/intermediate_data/fish_sp_matrix_stacked_array.rds")

