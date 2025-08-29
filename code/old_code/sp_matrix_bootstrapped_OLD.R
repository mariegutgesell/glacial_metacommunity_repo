##Species matrix -- generating arrays from bootstrapped samples 

library(tidyverse)
library(purrr)

##Read in benthic invert data
df <- read.csv("data/SEAK_Inverts_IndMeasure_AbunTakeMean_Repaired.csv") 

df$Date <- as.Date(df$Date, "%m/%d/%y")
df <- df %>% mutate(Month_Year = format(Date, "%Y-%m"))

##Calculate total biomass of each family for each date
df_summary <- df %>%
  select(StreamID, Month_Year, Date, Origin.x, Order.x, Family.x, Genus, Biomass_Ind) %>%
  filter(Origin.x == "Aquatic") %>%
  group_by(StreamID, Month_Year, Date, Origin.x, Order.x, Family.x) %>% ##family is lowest level ID for all taxa
  summarise_at(vars(Biomass_Ind), list(sp_total_biomass = sum)) %>%
  mutate(site_type = case_when(
    startsWith(StreamID, "Herbert") ~ "Glacial",
    startsWith(StreamID, "Steep") ~ "Snow",
    startsWith(StreamID, "Peter") ~ "Rain", 
    startsWith(StreamID, "Montana") ~ "Mixed",
  )) %>%
  ungroup() %>%
  filter(site_type != "Mixed")


# Ensure dataframe is in the correct format -- create array for variance partitioning function
df_summary <- df_summary %>%
  arrange(Family.x, Date, StreamID)

test_df <- df_summary %>%
  filter(site_type == "Rain") %>%
  select(site_type, Month_Year, Date, Family.x, sp_total_biomass)
head(test_df)


bootstrap_sample <- function(df) {
  df %>%
    sample_frac(replace = TRUE)
}

bootstrap_communities <- function(df, n_bootstrap = 50) {
  df %>%
    group_by(site_type, Date) %>%
    group_split() %>%
    map_dfr(~{
      map_dfr(1:n_bootstrap,function(x) {
        bootstrap_sample(.x) %>%
          mutate(Bootstrap_ID = x)
      })
    }) %>%
            
    ungroup()
}


test <- bootstrap_communities(test_df, n_bootstrap = 50)

ggplot(test, aes(x = Date, y = sp_total_biomass, group = Family.x, color = Family.x)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~Bootstrap_ID)

df_bootstrap <- bootstrap_communities(df_summary, n_bootstrap = 50)


##Create Arrays for each possible site combination, randomly taking a bootstrapped community 
sites <- unique(df_bootstrap$StreamID)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")
# Initialize an empty list to store the xtabs arrays
xtabs_list_collection <- list()

#Loop through each combination and create an xtabs array
set.seed(123)
for(n in 1:50) {
  xtabs_list_bs <- list()

  for (i in 1:nrow(unique_combinations)) {
    combo <- unique_combinations[i, ]
    selected_sites <- c(combo$site1, combo$site2, combo$site3)
  
   site_counts <- table(selected_sites)
  
    replicated_data <- data.frame()
  
  #For each site in the combination, replicate its data based on how many times it occurs
       for (site in names(site_counts)) {
           count <- site_counts[[site]]
          site_data <- df_bootstrap %>% filter(StreamID == site)
    
        selected_bootstrap_ids <- sample(unique(site_data$Bootstrap_ID), count, replace = FALSE)
    
    # Replicate the data according to how many times the site appears in the combination
       for (replicate_number in 1:count) {
      bootstrap_id <- selected_bootstrap_ids[replicate_number]
      
      replicated_site_data <- site_data %>%
        filter(Bootstrap_ID == bootstrap_id)
      replicated_site_data$Replicate_Type <- replicate_number
      
      replicated_data <- rbind(replicated_data, replicated_site_data) 
    }
    
  }
  replicated_data <- replicated_data %>%
    mutate(streamID_2 = paste(StreamID, Replicate_Type, sep = "_"))
  #Create an xtabs array for the replicated dataframe
  xtabs_array <- xtabs(sp_total_biomass ~ Family.x + Month_Year + streamID_2, data = replicated_data)
  
  xtabs_list_bs[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
  }
 xtabs_list_collection[[paste0("Iteration_", n)]] <- xtabs_list_bs 
}

