##Plotting Bootstrapped combinations of sites 

library(tidyverse)
library(ggplot2)


source("code/sp_matrix_bootstrapped_2.R")


##Expand final xtabs list into a df for plotting 
# Initialize an empty list to store the dataframes
list_of_dfs <- list()

# Iterate through each iteration (1 to 50)
for (iteration_num in seq_along(xtabs_list_collection)) {
  iteration <- xtabs_list_collection[[iteration_num]]  # Get the iteration
  
  # Iterate through each matrix in the iteration
  for (matrix_name in names(iteration)) {
    matrix_data <- iteration[[matrix_name]]  # Get the matrix
    
    # Convert the matrix to long-format dataframe
    long_df <- as.data.frame(matrix_data) %>%
      mutate(
        source_matrix = matrix_name,  # Add the matrix name
        iteration = iteration_num  # Add the iteration number
      )
    
    # Add the long-format dataframe to the list
    list_of_dfs[[length(list_of_dfs) + 1]] <- long_df
  }
}

# Bind all the dataframes together
combined_df <- bind_rows(list_of_dfs) %>%
  rename(biomass = "Freq")


##Calculate Total monthly biomass
combined_df_total <- combined_df %>%
  group_by(source_matrix, iteration, streamID_2, Month_Year) %>%
  summarise_at(vars(biomass), list(biomass_total = sum)) %>%
  group_by(source_matrix, streamID_2, Month_Year) %>%
  summarise_at(vars(biomass_total), list(biomass_total_mean = mean, biomass_total_sd = sd))


##All 3 Herbert Rivers
combined_df_total %>%
  filter(source_matrix == "Herbert River_Herbert River_Herbert River") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


combined_df_total %>%
  filter(source_matrix == "Herbert River_Herbert River_Peterson Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

combined_df_total %>%
  filter(source_matrix == "Herbert River_Herbert River_Steep Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


combined_df_total %>%
  filter(source_matrix == "Herbert River_Peterson Creek_Peterson Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

combined_df_total %>%
  filter(source_matrix == "Herbert River_Steep Creek_Steep Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


combined_df_total %>%
  filter(source_matrix == "Peterson Creek_Steep Creek_Steep Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_df_total %>%
  filter(source_matrix == "Peterson Creek_Peterson Creek_Steep Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


combined_df_total %>%
  filter(source_matrix == "Herbert River_Peterson Creek_Steep Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_df_total %>%
  filter(source_matrix == "Peterson Creek_Peterson Creek_Peterson Creek") %>%
  ggplot(aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(combined_df_total, aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~source_matrix)


###Total CV of all 3 communities 
total_all <- combined_df_total %>%
  group_by(source_matrix, Month_Year) %>%
  summarise_at(vars(biomass_total_mean), list(total = sum))

ggplot(total_all, aes(x = Month_Year, y = total, group = source_matrix, color = source_matrix)) +
  geom_point() +
  geom_line() +
 # geom_errorbar(aes(ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot() +
  geom_point(data = combined_df_total, aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_line(data = combined_df_total, aes(x = Month_Year, y = biomass_total_mean, group = streamID_2, color = streamID_2)) +
  geom_errorbar(data = combined_df_total, aes(x = Month_Year, ymin = biomass_total_mean - biomass_total_sd, ymax = biomass_total_mean + biomass_total_sd)) +
  geom_line(data= total_all, aes(x = Month_Year, y = total, group = source_matrix))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~source_matrix)
