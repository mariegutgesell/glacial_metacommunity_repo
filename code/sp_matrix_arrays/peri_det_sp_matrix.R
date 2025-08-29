##Create species matrix arrays for detritus and periphyton data

library(tidyverse)

peri_df <- read.csv("data/SEAK_Periphyton_Summary_Data.csv")
det_df <- read.csv("data/SEAK_Detritus_Summary_Data.csv")


##Add month_date to standardize to sampling month 
peri_df$Date <- as.Date(peri_df$Date, "%Y-%m-%d")
peri_df <- peri_df %>% mutate(Month_Year = format(Date, "%Y-%m")) %>%
  filter(Site != "Transitional") %>%
  select(Site:Month_Year) %>%
  mutate(Taxa = "Periphyton") %>%
  filter(Month_Year < "2019-01")
##for snow-fed, have two november sampling points, which one to use? 

det_df$Date <- as.Date(det_df$Date, "%Y-%m-%d")
det_df <- det_df %>% mutate(Month_Year = format(Date, "%Y-%m")) %>%
  filter(Site != "Transitional") %>%
  select(Site:Month_Year) %>%
  mutate(Taxa = "Detritus") %>%
  filter(Month_Year < "2019-01")


peri_matrix <- xtabs(meanAFDM ~ Taxa + Month_Year + Site, data = peri_df)
##does this make sense? don't have within stream different taxa... lets see if this works.. or put periphyton and detritus together in 1 dF?

det_matrix <- xtabs(meanAFDM_Area ~ Taxa + Month_Year + Site, data = det_df)
##are these two meanAFDM comparable? or is one corrected to area and the other not? 


##Create a list of arrays with all unique possible site combinations (not including mixed site)-- PERIPHYTON --------
sites <- unique(peri_df$Site)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")
# Initialize an empty list to store the xtabs arrays
peri_xtabs_list <- list()

#Loop through each combination and create an xtabs array
for (i in 1:nrow(unique_combinations)) {
  combo <- unique_combinations[i, ]
  selected_sites <- c(combo$site1, combo$site2, combo$site3)
  
  site_counts <- table(selected_sites)
  
  replicated_data <- data.frame()
  
  #For each site in the combination, replicate its data based on how many times it occurs
  for (site in names(site_counts)) {
    count <- site_counts[[site]]
    site_data <- peri_df[peri_df$Site == site, ]
    
    # Replicate the data according to how many times the site appears in the combination
    for (replicate_number in 1:count) {
      
      replicated_site_data <- site_data
      replicated_site_data$Replicate_Type <- replicate_number
      
      replicated_data <- rbind(replicated_data, replicated_site_data) 
    }
    
  }
  replicated_data <- replicated_data %>%
    mutate(streamID_2 = paste(Site, Replicate_Type, sep = "_"))
  #Create an xtabs array for the replicated dataframe
  xtabs_array <- xtabs(meanAFDM ~ Taxa + Month_Year + streamID_2, data = replicated_data)
  
  peri_xtabs_list[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
}

##need to fix so that all arrays have 3 levels 
saveRDS(peri_xtabs_list, "data/intermediate_data/peri_sp_matrix_stacked_array.rds")
##Create a list of arrays with all unique possible site combinations (not including mixed site)-- DETRITUS --------
sites <- unique(det_df$Site)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")
# Initialize an empty list to store the xtabs arrays
det_xtabs_list <- list()

#Loop through each combination and create an xtabs array
for (i in 1:nrow(unique_combinations)) {
  combo <- unique_combinations[i, ]
  selected_sites <- c(combo$site1, combo$site2, combo$site3)
  
  site_counts <- table(selected_sites)
  
  replicated_data <- data.frame()
  
  #For each site in the combination, replicate its data based on how many times it occurs
  for (site in names(site_counts)) {
    count <- site_counts[[site]]
    site_data <- det_df[det_df$Site == site, ]
    
    # Replicate the data according to how many times the site appears in the combination
    for (replicate_number in 1:count) {
      
      replicated_site_data <- site_data
      replicated_site_data$Replicate_Type <- replicate_number
      
      replicated_data <- rbind(replicated_data, replicated_site_data) 
    }
    
  }
  replicated_data <- replicated_data %>%
    mutate(streamID_2 = paste(Site, Replicate_Type, sep = "_"))
  #Create an xtabs array for the replicated dataframe
  xtabs_array <- xtabs(meanAFDM_Area ~ Taxa + Month_Year + streamID_2, data = replicated_data)
  
  det_xtabs_list[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
}

saveRDS(det_xtabs_list, "data/intermediate_data/det_sp_matrix_stacked_array.rds")

