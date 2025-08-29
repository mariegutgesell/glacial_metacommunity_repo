##plotting mean biomass combinations overtime

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)


inverts <- readRDS("data/intermediate_data/invert_sp_matrix_stacked_array_bootstrapped.rds")
peri <- readRDS("data/intermediate_data/peri_sp_matrix_stacked_array_bootstrapped.rds")

fish <- readRDS("data/intermediate_data/fish_sp_matrix_stacked_array_bootstrapped.rds")

str(inverts)

##Set up invert data -- randomly select 3 iterations, expand arrays into usable dfs
set.seed(123)
##randomly select 3 iterations, will use to plot 
select_iterations1 <- sample(1:50, 1)
selected_inverts <- inverts[select_iterations1]

#  Expand each of the 10 `xtabs` arrays into a dataframe
expand_xtabs_to_df <- function(xtabs_array, array_name, iteration_name) {
  df <- as.data.frame(as.table(xtabs_array)) %>%
    mutate(Combo = array_name,    # Retain the xtabs array name
           Iteration = iteration_name) # Track the original iteration
  return(df)
}

# loop through the selected elements and extract data
expanded_data <- lapply(seq_along(selected_inverts), function(i) {
  iteration_name <- names(selected_inverts)[i] # Get iteration name
  xtabs_list <- selected_inverts[[i]]         # Get the list of xtabs arrays
  
  # Expand each xtabs array into a dataframe
  df_list <- lapply(names(xtabs_list), function(array_name) {
    expand_xtabs_to_df(xtabs_list[[array_name]], array_name, iteration_name)
  })
  
  return(bind_rows(df_list)) # Combine all expanded dataframes for this iteration
})

# combine the results into one dataframe
final_df <- bind_rows(expanded_data)

inverts_df <- final_df %>%
  filter(Combo %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  mutate(combo_type = case_when(
    startsWith(Combo, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combo, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combo, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
  )) %>%
  mutate(taxa = "inverts") %>%
  group_by(taxa, Combo, combo_type, streamID_2, Iteration, Month_Year) %>%
  summarise_at(vars(Freq), list(iteration_total_biomass = sum)) %>%
  separate(streamID_2, into = c("stream_type", "rep"), sep = "_", remove =FALSE)
  
inverts_df_total_biomass <- inverts_df %>%  
  group_by(taxa, Combo, combo_type, Month_Year) %>%
  summarise_at(vars(iteration_total_biomass), list(total_biomass = sum)) %>%
  rename(Combination = "Combo")

inverts_df_total_biomass$combo_type <- ordered(inverts_df_total_biomass$combo_type,
                                            levels = c("Heterogenous",  "2:1","Homogenous" )) 

inverts_df$combo_type <- ordered(inverts_df$combo_type,
                              levels = c("Heterogenous",  "2:1","Homogenous" )) 



inverts_df$month <- months(as.Date(paste("01", inverts_df$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)
inverts_df_total_biomass$month <- months(as.Date(paste("01", inverts_df_total_biomass$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)


##Set up peri data -- randomly select 3 iterations, expand arrays into usable dfs
set.seed(123)
##randomly select 3 iterations, will use to plot 
select_iterations1 <- sample(1:50, 1)
selected_peri <- peri[select_iterations1]


# loop through the selected elements and extract data
expanded_data <- lapply(seq_along(selected_peri), function(i) {
  iteration_name <- names(selected_peri)[i] # Get iteration name
  xtabs_list <- selected_peri[[i]]         # Get the list of xtabs arrays
  
  # Expand each xtabs array into a dataframe
  df_list <- lapply(names(xtabs_list), function(array_name) {
    expand_xtabs_to_df(xtabs_list[[array_name]], array_name, iteration_name)
  })
  
  return(bind_rows(df_list)) # Combine all expanded dataframes for this iteration
})

#combine the results into one dataframe
final_df <- bind_rows(expanded_data)

peri_df <- final_df %>%
  filter(Combo %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  mutate(combo_type = case_when(
    startsWith(Combo, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combo, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combo, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
  )) %>%
  mutate(taxa = "Periphyton") %>%
  group_by(taxa, Combo, combo_type, streamID_2, Iteration, Month_Year) %>%
  summarise_at(vars(Freq), list(iteration_total_biomass = sum)) %>%
  separate(streamID_2, into = c("stream_type", "rep"), sep = "_", remove =FALSE)

peri_df_total_biomass <- peri_df %>%  
  group_by(taxa, Combo, combo_type, Month_Year) %>%
  summarise_at(vars(iteration_total_biomass), list(total_biomass = sum)) %>%
  rename(Combination = "Combo")


peri_df_total_biomass$combo_type <- ordered(peri_df_total_biomass$combo_type,
                                            levels = c("Heterogenous",  "2:1","Homogenous" )) 

peri_df$combo_type <- ordered(peri_df$combo_type,
                              levels = c("Heterogenous",  "2:1","Homogenous" )) 



peri_df$month <- months(as.Date(paste("01", peri_df$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)
peri_df_total_biomass$month <- months(as.Date(paste("01", peri_df_total_biomass$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)

##Set up fish data -- randomly select 3 iterations, expand arrays into usable dfs
set.seed(123)
##randomly select 3 iterations, will use to plot 
select_iterations1 <- sample(1:50, 1)
selected_fish <- fish[select_iterations1]


# loop through the selected elements and extract data
expanded_data <- lapply(seq_along(selected_fish), function(i) {
  iteration_name <- names(selected_fish)[i] # Get iteration name
  xtabs_list <- selected_fish[[i]]         # Get the list of xtabs arrays
  
  # Expand each xtabs array into a dataframe
  df_list <- lapply(names(xtabs_list), function(array_name) {
    expand_xtabs_to_df(xtabs_list[[array_name]], array_name, iteration_name)
  })
  
  return(bind_rows(df_list)) # Combine all expanded dataframes for this iteration
})

# combine the results into one dataframe
final_df <- bind_rows(expanded_data)

fish_df <- final_df %>%
  filter(Combo %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  mutate(combo_type = case_when(
    startsWith(Combo, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combo, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combo, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
  )) %>% 
  mutate(taxa = "fish") %>%
  group_by(taxa, Combo, combo_type, streamID_2, Iteration, Month_Year) %>%
  summarise_at(vars(Freq), list(iteration_total_biomass = sum)) %>%
  separate(streamID_2, into = c("stream_type", "rep"), sep = "_", remove =FALSE)

fish_df_total_biomass <- fish_df %>%  
  group_by(taxa, Combo, combo_type, Month_Year) %>%
  summarise_at(vars(iteration_total_biomass), list(total_biomass = sum)) %>%
  rename(Combination = "Combo")

fish_df_total_biomass$combo_type <- ordered(fish_df_total_biomass$combo_type,
                                       levels = c("Heterogenous",  "2:1","Homogenous" )) 

fish_df$combo_type <- ordered(fish_df$combo_type,
                                            levels = c("Heterogenous",  "2:1","Homogenous" )) 



fish_df$month <- months(as.Date(paste("01", fish_df$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)
fish_df_total_biomass$month <- months(as.Date(paste("01", fish_df_total_biomass$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)




##Plotting out biomass timeseries
##Fish 

fish_plot <- ggplot()+
  geom_line(data = fish_df, aes(x = Month_Year, y = iteration_total_biomass, group = streamID_2, color = site_type, linetype = stream_type),color = "darksalmon", alpha = 0.6) +
  geom_line(data = fish_df_total_biomass , aes(x = Month_Year, y = total_biomass, group = combo_type, color = combo_type),color = "darksalmon", linewidth=1) +
  theme_classic() +
  xlab("Date") +
  ylab("Total Biomass (g)") +
 # ylim(0,2700)+
  theme(axis.title.x = element_blank(), text = element_text(family = "Times New Roman", size = 14), legend.position  = "none", strip.text = element_blank()) +
  facet_wrap(~combo_type) +
  scale_x_discrete(labels = fish_df$month)

fish_plot



inverts_plot <- ggplot()+
  geom_line(data = inverts_df, aes(x = Month_Year, y = iteration_total_biomass, group = streamID_2, color = site_type, linetype = stream_type),color = "deepskyblue", alpha = 0.6) +
  geom_line(data = inverts_df_total_biomass , aes(x = Month_Year, y = total_biomass, group = combo_type, color = combo_type),color = "deepskyblue", linewidth=1) +
  theme_classic() +
  xlab("Date") +
  ylab("Total Biomass (g/m2)") +
  # ylim(0,2700)+
  theme(axis.title.x = element_blank(), text = element_text(family = "Times New Roman", size = 14), legend.position  = "none", strip.text = element_blank()) +
  facet_wrap(~combo_type) +
  scale_x_discrete(labels = inverts_df$month)

inverts_plot



peri_plot <- ggplot()+
  geom_line(data = peri_df, aes(x = Month_Year, y = iteration_total_biomass, group = streamID_2, color = site_type, linetype = stream_type),color = "darkgreen", alpha = 0.6) +
  geom_line(data = peri_df_total_biomass , aes(x = Month_Year, y = total_biomass, group = combo_type, color = combo_type),color = "darkgreen", linewidth=1) +
  theme_classic() +
  xlab("Date") +
  ylab("Total Biomass (g/m2)") +
  # ylim(0,2700)+
  theme(axis.title.x = element_blank(), text = element_text(family = "Times New Roman", size = 14), legend.position  = "none", strip.text = element_blank()) +
  facet_wrap(~combo_type) +
  scale_x_discrete(labels = peri_df$month)

peri_plot


all_plot<- ggarrange(peri_plot, inverts_plot, fish_plot,
                     
                     ncol = 1, nrow = 3, labels = c("a)", "b)", "c)"), font.label = list(colour = "black", size = 14, family = "Times New Roman"))
all_plot

