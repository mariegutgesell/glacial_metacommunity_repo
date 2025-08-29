##Code to create site species matrix of benthic invertebrate biomass data -- no bootstrapping
##M.Gutgesell
##Nov 13, 2024

library(tidyverse)
library(lubridate)

##Read in benthic invert data
df <- read.csv("data/SEAK_Inverts_IndMeasure_AbunTakeMean_Repaired.csv") 
  
##note: Biomass = a*BodyLength^b - calculated at individual level 
df$Date <- as.Date(df$Date, "%m/%d/%y")
df <- df %>% mutate(Month_Year = format(Date, "%Y-%m"))

df_na <- df %>%
  filter(is.na(Biomass_Ind)) 


test <- df %>%
  filter(StreamID == "Steep Creek", Month_Year == "2018-04") %>%
  group_by(Sample_Number, Order.x, Family.x) %>%
  count()



##Calculate total biomass of each family for each date
df_summary <- df %>%
  select(StreamID, Month_Year, Date, Origin.x, Order.x, Family.x, Genus, Biomass_Ind) %>%
  filter(Origin.x == "Aquatic") %>%
  group_by(StreamID, Month_Year, Date, Origin.x, Order.x, Family.x) %>% ##family is lowest level ID for all taxa
  summarise_at(vars(Biomass_Ind), list(sp_total_biomass = sum)) %>%
  mutate(site_type = case_when(
    startsWith(StreamID, "Herbert") ~ "Glacier-fed",
    startsWith(StreamID, "Steep") ~ "Snow-fed",
    startsWith(StreamID, "Peter") ~ "Rain-fed", 
    startsWith(StreamID, "Montana") ~ "Mixed",
  )) %>%
  ungroup() %>%
  filter(Month_Year < "2019-01")

df_na <- df_summary %>%
  filter(is.na(sp_total_biomass))
# Ensure dataframe is in the correct format -- create array for variance partitioning function
df_summary <- df_summary %>%
  arrange(Family.x, Month_Year, site_type)
# Create the array
sp_matrix_all <- xtabs(sp_total_biomass ~ Family.x + Month_Year + site_type, data = df_summary)


##Create array w/ mixed site removed
df_summary_2 <- df_summary %>%
  filter(site_type != "Mixed")
df_summary_2 <- df_summary_2 %>%
  arrange(Family.x, Date, site_type)
sp_matrix_2 <- xtabs(sp_total_biomass ~ Family.x + Month_Year + site_type, data = df_summary_2)


saveRDS(sp_matrix_2, "data/intermediate_data/invert_sp_matrix.rds")
#sp_matrix <- readRDS("data/intermediate_data/sp_matrix.rds")
str(sp_matrix_2)



##Create a list of arrays with all unique possible site combinations (not including mixed site) --------
sites <- unique(df_summary_2$site_type)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")
# Initialize an empty list to store the xtabs arrays
xtabs_list <- list()

#Loop through each combination and create an xtabs array
for (i in 1:nrow(unique_combinations)) {
  combo <- unique_combinations[i, ]
  selected_sites <- c(combo$site1, combo$site2, combo$site3)

  site_counts <- table(selected_sites)
  
  replicated_data <- data.frame()
  
  #For each site in the combination, replicate its data based on how many times it occurs
  for (site in names(site_counts)) {
    count <- site_counts[[site]]
    site_data <- df_summary_2[df_summary_2$site_type == site, ]
    
    # Replicate the data according to how many times the site appears in the combination
    for (replicate_number in 1:count) {

      replicated_site_data <- site_data
      replicated_site_data$Replicate_Type <- replicate_number

      replicated_data <- rbind(replicated_data, replicated_site_data) 
    }

  }
  replicated_data <- replicated_data %>%
    mutate(streamID_2 = paste(site_type, Replicate_Type, sep = "_"))
  #Create an xtabs array for the replicated dataframe
  xtabs_array <- xtabs(sp_total_biomass ~ Family.x + Month_Year + streamID_2, data = replicated_data)
  
  xtabs_list[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
}

saveRDS(xtabs_list, "data/intermediate_data/invert_sp_matrix_stacked_array.rds")


