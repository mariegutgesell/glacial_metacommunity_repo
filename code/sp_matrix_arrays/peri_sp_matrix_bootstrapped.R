##Create species matrix arrays for periphyton data -- attempting bootstrap

library(tidyverse)

peri_df_summary <- read.csv("data/SEAK_Periphyton_Summary_Data.csv")

peri_df <- read.csv("data/SEAK_Periphyton.csv")


##Clean up dateframe and insert missing dates associated with a specific sample code, if present 
##Add month_date to standardize to sampling month 
peri_df$Date <- as.Date(peri_df$Date, "%Y-%m-%d")
peri_df_clean <- peri_df %>% 
  mutate(base_sample_code = str_extract(SampCode, "^(.*?)(?=Peri)")) %>% 
  select(Site, Date, SampleID, base_sample_code, SampCode, AFDM, Area_cm2:SubsampleArea_cm2) %>%
  group_by(base_sample_code) %>%
  fill(Date, .direction = "downup") %>% ##fill in NA dates with associated date from that same sample (some replicates from a given date had missing sample date)
  mutate(Month_Year = format(Date, "%Y-%m")) %>%
  filter(Site != "Montana") %>%
  select(Site:Month_Year) %>%
  mutate(Taxa = "Periphyton") %>%
  filter(Month_Year < "2019-01") %>%
  mutate(site_type = case_when(
    startsWith(Site, "Herbert") ~ "Glacier-fed",
    startsWith(Site, "Steep") ~ "Snow-fed",
    startsWith(Site, "Peter") ~ "Rain-fed", 
  ))


duplicates <- peri_df_clean[duplicated(peri_df_clean) | duplicated(peri_df_clean, fromLast = TRUE), ]
##there are two Peri_04 for april Herbert sample that are exactly the same, for now will remove that duplicate... 

peri_df_clean <- peri_df_clean %>%
  distinct()

##Scale AFDM to m2 
df_scaled <- peri_df_clean %>%
  select(Site, site_type, Taxa, Date, Month_Year, SampleID, base_sample_code, SampCode, AFDM, Area_cm2:SubsampleArea_cm2) %>%
  mutate(AFDM_mg_per_cm2 = AFDM/SubsampleArea_cm2) %>%
  mutate(AFDM_mg_per_m2 = AFDM_mg_per_cm2*1000) %>%
  mutate(AFDM_g_per_m2 = AFDM_mg_per_m2/1000)


test <- df_scaled %>%
  filter(Site == "Herbert")

##Create bootstrap dataframes, selecting 5 peri samples randomly (with replacement, so always have 5, but different combinations of 5)
bootstrap_sample <- function(df) {
  selected_samples <- sample(unique(df$SampleID), 5, replace = TRUE)
  sampled_df <- do.call(
    rbind,
    lapply(selected_samples, function(sample_num) {
      df %>% filter(SampleID == sample_num)
    })
  )
  
  return(sampled_df)
}

bootstrap_samples <- function(df, n_bootstrap = 50) {
  df %>%
    group_by(site_type, Month_Year) %>%
    group_split() %>%
    map_dfr(~{
      map_dfr(1:n_bootstrap,function(x) {
        bootstrap_sample(.x) %>%
          mutate(Bootstrap_ID = x)
      })
    }) %>%
    
    ungroup()
}


df_bootstrap <- bootstrap_samples(df_scaled, n_bootstrap = 50)



test <- df_bootstrap %>%
  filter(Bootstrap_ID ==1)

##Calculate meanAFDM for each bootstrap ID 
df_bootstrap_mean <- df_bootstrap %>%
  group_by(Taxa, Site, site_type, Date, Month_Year, Bootstrap_ID) %>%
  summarise_at(vars(AFDM_g_per_m2), list(meanAFDM_g_per_m2 = mean))




##Create arrays for each bootstrapped community for each combination of sites
sites <- unique(df_bootstrap_mean$site_type)
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
      site_data <- df_bootstrap_mean %>% filter(site_type == site)
      
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
      xtabs_array <- xtabs(meanAFDM_g_per_m2 ~ Taxa + Month_Year + streamID_2, data = replicated_data)
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



saveRDS(xtabs_list_collection, "data/intermediate_data/peri_sp_matrix_stacked_array_bootstrapped.rds")

