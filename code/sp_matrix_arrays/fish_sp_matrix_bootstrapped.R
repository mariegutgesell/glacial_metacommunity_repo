##Create bootstrapped species matrix arrays for fish 
##M.Gutgesell
##Feb 18, 2025

library(tidyverse)
library(readxl)

##Import fish data - already standardized to area/effort, use mean and 95% CI of biomass estimate to create normal distribution and then sample from that to bootstrap
fish_df <- read_excel("data/SEAK_Fish_Biomass_Supplement.xlsx")
fish_df$Date <- as.Date(fish_df$Date, "%m-%d-%y")

fish_df <- fish_df %>%   mutate(Month_Year = format(Date, "%Y-%m")) %>%
  filter(Stream != "Transitional") %>%
  mutate(species2 = str_extract(Species, "^(.*?)(?=Juvenile)")) %>% 
  mutate(species3 = str_extract(Species, "^(.*?)(?=Fry)")) %>% 
  mutate(species2 = coalesce(species2, species3)) %>%
  mutate(species2 = coalesce(species2, Species)) %>%
  separate(`95% CI`, into = c("lower_95_CI", "upper_95_CI"), sep = "-")

fish_df$lower_95_CI <- as.numeric(fish_df$lower_95_CI)
fish_df$upper_95_CI <- as.numeric(fish_df$upper_95_CI)


##For each species and each sampling point per site, create a distribution using the mean and 95% CI 
##first calcuate std dev for each
fish_df<- fish_df %>%
  mutate(sd_1 = (upper_95_CI - Biomass)/1.96) %>% ##this assumes mean is in middle... 
  mutate(sd_2 = (upper_95_CI - lower_95_CI) / (2 * 1.96)) 

head(fish_df)


n_samples <- 1000

fish_df <- fish_df %>%
  mutate(distribution = map2(Biomass, sd_2, ~ rnorm(n_samples, mean = .x, sd = .y))) %>%
  filter(Month_Year < "2019-01")

##Visualizing to see if have distributions right
df_long <- fish_df%>%
  unnest_longer(distribution)

df_long %>%
  filter(Stream == "Snow-fed", species2 == "DollyVarden", Month_Year == "2018-08") %>% 
  ggplot(aes(x = distribution, fill = as.factor(Month_Year))) +
  geom_density() +
  facet_wrap(~Month_Year, scale = "free") +
  theme_classic()

df_long %>%
  filter(Stream == "Rain-fed", species2 == "Coastal Trout") %>% 
  ggplot(aes(x = distribution, fill = as.factor(Month_Year))) +
  geom_density() +
  facet_wrap(~Month_Year, scale = "free") +
  theme_classic()

##separate fry/juvenile or not? ... ? - no, only glacial stream the juvenile and fry is part of name, no site/date where both occur, use name without juvenile so comparable across sites 





fish_matrix <- xtabs(Biomass ~ species2 + Month_Year + Stream, data = fish_df)


##need to bootstrap data...so for each replicate in the 3 streams randomly takes a different biomass, and also creates 50 different possible stacked arrays... 



##Create a list of arrays with all unique possible site combinations (not including mixed site)-- Fish --------
set.seed(123)
sites <- unique(fish_df$Stream)
site_combinations <- expand.grid(site1 = sites, site2 = sites, site3 = sites)
site_combinations_sorted <- t(apply(site_combinations, 1, sort))
unique_combinations <- unique(as.data.frame(site_combinations_sorted))
colnames(unique_combinations) <- c("site1", "site2", "site3")
# Initialize an empty list to store the xtabs arrays
fish_xtabs_list_collection <- list()

#num_bootstrap <- 50

#Loop through each combination and create an xtabs array
# Loop to create 50 arrays using bootstrapped communities for each site combination
for (n in 1:50) {
  xtabs_list_bs <- list()

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
      boot_data <- site_data %>%
       # group_by(Stream, Species, Month_Year) %>%
        mutate(distribution = map2(Biomass, sd_2, ~ rnorm(n_samples, mean = .x, sd = .y))) %>%
        mutate(bootstrapped_biomass = map_dbl(distribution, ~ sample(.x, 1))) %>%  # Select a single value
        select(-distribution)
     
      boot_data$Stream <- site
      boot_data$Replicate_Type <- replicate_number
      replicated_data <- bind_rows(replicated_data, boot_data)
    }
  }
  
  replicated_data <- replicated_data %>%
    mutate(streamID_2 = paste(Stream, Replicate_Type, sep = "_"))
  
  if (nrow(replicated_data) > 0) {
    xtabs_array <- xtabs(bootstrapped_biomass ~ species2 + Month_Year + streamID_2, data = replicated_data)
    xtabs_list_bs[[paste(combo$site1, combo$site2, combo$site3, sep = "_")]] <- xtabs_array
    }
  }
  
  # Add the list of xtabs arrays for this iteration to the collection
  fish_xtabs_list_collection[[paste0("Iteration_", n)]] <- xtabs_list_bs 
}

negative_vals <- xtabs_array[xtabs_array < 0]


saveRDS(fish_xtabs_list_collection, "data/intermediate_data/fish_sp_matrix_stacked_array_bootstrapped.rds")



##testing out lognormal distribution
# Given mean and sd on original scale
m <- 178.954706
s <- 105.905762

# Compute lognormal parameters
sigma_log2 <- log(1 + (s/m)^2)
sigma_log <- sqrt(sigma_log2)
mu_log <- log(m) - 0.5 * sigma_log2

# Generate values
set.seed(123)
samples <- rlnorm(1000, meanlog = mu_log, sdlog = sigma_log)

# Optional: check the mean and sd of simulated values
mean(samples)
sd(samples)

samples_df <- as.data.frame(samples)

samples_df %>%
#  filter(Stream == "Snow-fed", species2 == "DollyVarden") %>% 
  ggplot(aes(x = samples)) +
  geom_density() +
#  facet_wrap(~Month_Year, scale = "free") +
  theme_classic()
