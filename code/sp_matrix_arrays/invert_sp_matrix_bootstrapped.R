##Species matrix -- generating arrays from bootstrapped samples -- bootstrapped based on samples -- RIGHT CODE

library(tidyverse)
library(purrr)

##Read in benthic invert data
df <- read.csv("data/SEAK_Inverts_IndMeasure_AbunTakeMean_Repaired.csv") 

df$Date <- as.Date(df$Date, "%m/%d/%y")
df <- df %>% mutate(Month_Year = format(Date, "%Y-%m")) %>%
  mutate(site_type = case_when(
    startsWith(StreamID, "Herbert") ~ "Glacier-fed",
    startsWith(StreamID, "Steep") ~ "Snow-fed",
    startsWith(StreamID, "Peter") ~ "Rain-fed", 
    startsWith(StreamID, "Montana") ~ "Mixed",
  )) %>%
  ungroup() %>%
  filter(site_type != "Mixed") %>%
  filter(Month_Year < "2019-01")


df$Split <- gsub("\\*\\*", "", df$Split) 
##what does the 1** mean for herbert surber_03 2018-11-12?  - have removed the ** for now 
df$Split <- as.numeric(df$Split)

##NEED TO STANDARDIZE BIOMASS TO M2
##I think this data is just individual level biomass estimates for everything in sample 
##need to standardize to m2 based on subsample volume and area of benthic habitat sampled -- need this information 
##note: Biomass = a*BodyLength^b - calculated at individual level - but i think still need to scale this to m2 

##Standardizing biomass
df_summary <- df %>%
  select(site_type, StreamID, Month_Year, Date, Sample_Number, Split, Origin.x, Order.x, Family.x, Genus, Taxon, Biomass_Ind, Large_Rare) %>%
  group_by(site_type, StreamID, Month_Year, Date, Sample_Number, Split, Origin.x, Order.x, Family.x, Taxon, Large_Rare) %>% ##family is lowest level ID for all taxa
  summarise(
    sp_biomass_subsample = sum(Biomass_Ind, na.rm = TRUE),
    n_individuals = n()
  ) %>%
  mutate(subsample_multiplier = 1/Split) %>%
  mutate(sp_total_biomass = sp_biomass_subsample*subsample_multiplier,
         n_individuals_total = n_individuals*subsample_multiplier) %>%
  group_by(site_type, StreamID, Month_Year, Date, Sample_Number, Split, Origin.x, Order.x, Family.x, Taxon) %>%
  summarise(
    sp_total_biomass_mg = sum(sp_total_biomass, na.rm = TRUE),
    n_individuals_total = sum(n_individuals_total)
  ) %>%
  mutate(sp_total_biomass_g_m2 = (sp_total_biomass_mg/1000)*1.33)




##Create bootstrap dataframes, selecting 5 surber samples randomly 
bootstrap_sample <- function(df) {
  selected_samples <- sample(unique(df$Sample_Number), 5, replace = TRUE)
  sampled_df <- do.call(
    rbind,
    lapply(selected_samples, function(sample_num) {
      df %>% filter(Sample_Number == sample_num)
    })
  )
  
  return(sampled_df)
}

bootstrap_samples <- function(df, n_bootstrap = 50) {
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


df_bootstrap <- bootstrap_samples(df_summary, n_bootstrap = 50)

test <- df_bootstrap %>%
  filter(Bootstrap_ID ==1)

##Create arrays for each bootstrapped community for each combination of sites
sites <- unique(df_bootstrap$site_type)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")

# Initialize an empty list to store the xtabs arrays
xtabs_list_collection <- list()

#Loop through each combination and create an xtabs array
set.seed(123)  # Set seed for reproducibility

# Loop to create 50 arrays using bootstrapped communities for each site combination
for (n in 1:50) {
  xtabs_list_bs <- list()
  
  # Loop through each unique combination of sites
  for (i in 1:nrow(unique_combinations)) {
    combo <- unique_combinations[i, ]
    selected_sites <- c(combo$site1, combo$site2, combo$site3)
    
    # Count how many times each site appears in the combination
    site_counts <- table(selected_sites)
    
    # Initialize an empty dataframe to store the replicated data
    replicated_data <- data.frame()
    
    # Loop through each site in the combination
    for (site in names(site_counts)) {
      count <- site_counts[[site]]
      
      # Filter the data for the current site
      site_data <- df_bootstrap %>% filter(site_type == site)
      
      # Select unique Bootstrap_IDs for each repetition of the site
      if (nrow(site_data) > 0 && length(unique(site_data$Bootstrap_ID)) >= count) {
        selected_bootstrap_ids <- sample(unique(site_data$Bootstrap_ID), count, replace = FALSE)
      } else {
        selected_bootstrap_ids <- sample(unique(site_data$Bootstrap_ID), count, replace = TRUE)
      }
      
      # Replicate the data for the selected Bootstrap_IDs
      for (replicate_number in 1:count) {
        bootstrap_id <- selected_bootstrap_ids[replicate_number]
        
        replicated_site_data <- site_data %>%
          filter(Bootstrap_ID == bootstrap_id) %>%
          mutate(Replicate_Type = replicate_number)
        
        replicated_data <- rbind(replicated_data, replicated_site_data)
      }
    }
    
    # Create a unique identifier for the combination of site and replicate
    replicated_data <- replicated_data %>%
      mutate(streamID_2 = paste(site_type, Replicate_Type, sep = "_"))
    
    # Create an xtabs array for the replicated dataframe
    if (nrow(replicated_data) > 0) {
      xtabs_array <- xtabs(sp_total_biomass_g_m2 ~ Family.x + Month_Year + streamID_2, data = replicated_data)
      xtabs_list_bs[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
    }
  }
  
  # Add the list of xtabs arrays for this iteration to the collection
  xtabs_list_collection[[paste0("Iteration_", n)]] <- xtabs_list_bs 
}


##Check that there are non-empty xtabs
non_empty_xtabs <- lapply(xtabs_list_collection, function(iteration_list) {
  iteration_list[!sapply(iteration_list, is.null)]
})


non_empty_xtabs <- non_empty_xtabs[sapply(non_empty_xtabs, length) > 0]



saveRDS(xtabs_list_collection, "data/intermediate_data/invert_sp_matrix_stacked_array_bootstrapped.rds")


