##Partition Metacommunity Variation


source("code/sp_matrix_arrays/invert_sp_matrix_bootstrapped.R")

rm(list = ls()[!ls() %in% c("xtabs_list_collection")])

##testing out variance partitioning code (from Wang et al., 2019, Ecography, 42: 1-12)

var.partition <- function(metacomm_tsdata){
  ## The function "var.partition" performs the partitioning of variability
  ## across hierarchical levesl within a metacommunity.
  ## The input array "metacomm_tsdata" is an N*T*M array. The first dimension represents N species,
  ## the second represents time-series observations of length T, and the third represents M local communities.
  ## The output includes four variability and four synchrony metrics as defined in the main text.
  ## Note that, to be able to handle large metacommunities, this code has avoided calculating all covariance.
  ts_metacom <- apply(metacomm_tsdata,2,sum)
  ts_patch <- apply(metacomm_tsdata,c(2,3),sum)
  ts_species <- apply(metacomm_tsdata,c(1,2),sum)
  sd_metacom <- sd(ts_metacom)
  sd_patch_k <- apply(ts_patch,2,sd)
  sd_species_i <- apply(ts_species,1,sd)
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd)
  mean_metacom <- mean(ts_metacom)
  CV_S_L <- sum(sd_species_patch_ik)/mean_metacom
  CV_C_L <- sum(sd_patch_k)/mean_metacom
  CV_S_R <- sum(sd_species_i)/mean_metacom
  CV_C_R <- sd_metacom/mean_metacom
  phi_S_L2R <- CV_S_R/CV_S_L
  phi_C_L2R <- CV_C_R/CV_C_L
  phi_S2C_L <- CV_C_L/CV_S_L
  phi_S2C_R <- CV_C_R/CV_S_R
  partition_3level <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R,
                        phi_S_L2R=phi_S_L2R, phi_C_L2R=phi_C_L2R, phi_S2C_L=phi_S2C_L,
                        phi_S2C_R=phi_S2C_R,
                        mean_metacom = mean_metacom, sd_metacom = sd_metacom)
  return(partition_3level)
}





##calculate CV and synchrony for each of the site combinations

##trying to get partitoning results for each of the 50 different arrays
partition_results <- list()

for(n in 1:50) {
  current_xtabs_list <- xtabs_list_collection[[paste0("Iteration_", n)]]
  
  iteration_results <- list()
    for(combo in names(current_xtabs_list)) {
      df <- current_xtabs_list[[combo]]
      partition_result <- var.partition(df)
      iteration_results[[combo]] <- partition_result
    }
  partition_results[[paste0("Iteration_", n)]] <- iteration_results
  
  
}


##summarise mean and standard error 
# Define a helper function to calculate standard error
standard_error <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

# Initialize a list to store the results for each combination and each variable
partition_grouped_summary <- list()

# Extract all unique combinations of site pairs from the first iteration
site_combinations <- names(partition_results[["Iteration_1"]])

# Extract the variable names from the first partition result
partition_variables <- names(partition_results[["Iteration_1"]][[1]])

# Loop through each site combination
for (combo in site_combinations) {
  
  # Initialize a list to store the mean and standard error for each variable
  combo_summary <- list()
  
  # Loop through each variable and calculate the mean and standard error across the 50 iterations for this combo
  for (var in partition_variables) {
    
    # Initialize a vector to store the values of this variable across the 50 iterations for the current combo
    var_values <- numeric(50)
    
    # Loop through each iteration to extract the value for the current variable for this combination
    for (n in 1:50) {
      # Extract the partition results for this iteration
      iteration_results <- partition_results[[paste0("Iteration_", n)]]
      
      # Get the value for the current combo and variable in this iteration
      var_values[n] <- iteration_results[[combo]][[var]]
    }
    
    # Calculate the mean and standard error for this variable across the 50 iterations
    mean_var <- mean(var_values, na.rm = TRUE)
    se_var <- standard_error(var_values)
    
    # Store the results for this variable in the combo_summary
    combo_summary[[var]] <- list(mean = mean_var, standard_error = se_var)
  }
  
  # Store the summary for this combo in the partition_grouped_summary list
  partition_grouped_summary[[combo]] <- combo_summary
}

# Convert the summary results to a data frame for better presentation
partition_grouped_summary_df <- do.call(rbind, lapply(partition_grouped_summary, function(combo) {
  combo_df <- do.call(rbind, lapply(combo, as.data.frame))
  combo_df$Variable <- rownames(combo_df)
  return(combo_df)
}))

# Add a column for the combination names
partition_grouped_summary_df$Combo <- rep(site_combinations, each = length(partition_variables))

# Rearrange columns for better readability
partition_grouped_summary_df <- partition_grouped_summary_df[, c("Combo", "Variable", "mean", "standard_error")]


##change to wide
summary_df <- partition_grouped_summary_df %>%
  pivot_wider(names_from = Variable, values_from = c(mean, standard_error), names_glue = "{Variable}_{.value}") %>%
  mutate(combo_type = case_when(
    startsWith(Combo, "Glacier-fed_Glacier-fed_Glacier-fed") ~ "Homogenous",
    startsWith(Combo, "Glacier-fed_Glacier-fed_Rain-fed") ~ "2:1",
    startsWith(Combo, "Glacier-fed_Glacier-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Glacier-fed_Rain-fed_Rain-fed") ~ "2:1",
    startsWith(Combo, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combo, "Glacier-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combo, "Rain-fed_Rain-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Snow-fed_Snow-fed_Snow-fed") ~ "Homogenous",
  ))
  
  

##stabilization from community to metacommunity CV
cv_diff <- summary_df %>%
  mutate(comm_cv_diff = CV_C_L_mean - CV_C_R_mean,
         pop_cv_diff = CV_S_L_mean - CV_S_R_mean)

##see greatest difference in CV from community to metacommunity and population to metapopulation for most heterogeneous combo
##also see greatest asynchrony 

##Plotting out some CV/synchrony results
cv_diff$combo_type <- ordered(cv_diff$combo_type, levels = c("Homogenous", "2:1", "Heterogenous"))



##Metacommunity Stability 
ggplot(cv_diff, aes(x = combo_type, y = CV_C_R_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Metacommunity Variability (CV-C,R) ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##Local Community Stability 
ggplot(cv_diff, aes(x = combo_type, y = CV_C_L_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = CV_C_L_mean - CV_C_L_standard_error, ymax = CV_C_L_mean + CV_C_L_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Local Community Variability (CV-C,L) ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##CV dampening
ggplot(cv_diff, aes(x = combo_type, y = comm_cv_diff, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
#  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Community to Metacommunity CV Dampening ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##Local population Variability
##weighted average of local population CV across communities
ggplot(cv_diff, aes(x = combo_type, y = CV_S_L_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = CV_S_L_mean - CV_S_L_standard_error, ymax = CV_S_L_mean + CV_S_L_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Local Population Variability")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##Metapopulation Variability
##weighted average CV of total metapopulation biomass
ggplot(cv_diff, aes(x = combo_type, y = CV_S_R_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = CV_S_R_mean - CV_S_R_standard_error, ymax = CV_S_R_mean + CV_S_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Metapopulation Variability")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##Asynchrony -- community asynchrony local to regional -- driven by spatial community dissimilarity 
##Is the spatial synchrony of total community biomass among local patches
ggplot(cv_diff, aes(x = combo_type, y = phi_C_L2R_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = phi_C_L2R_mean - phi_C_L2R_standard_error, ymax = phi_C_L2R_mean + phi_C_L2R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Community Level Spatial Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##Asynchrony -- local scale species synchrony
##
ggplot(cv_diff, aes(x = combo_type, y = phi_S2C_L_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = phi_S2C_L_mean - phi_S2C_L_standard_error, ymax = phi_S2C_L_mean + phi_S2C_L_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Local Scale Species Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##Asynchrony -- species level spatial synchrony 
##
ggplot(cv_diff, aes(x = combo_type, y = phi_S_L2R_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = phi_S_L2R_mean - phi_S_L2R_standard_error, ymax = phi_S_L2R_mean + phi_C_L2R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Species Level Spatial Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##Asynchrony -- regional scale species synchrony
##
ggplot(cv_diff, aes(x = combo_type, y = phi_S2C_R_mean, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = phi_S2C_R_mean - phi_S2C_R_standard_error, ymax = phi_S2C_R_mean + phi_S2C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Regional Scale Species Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##Population level CV dampening -- from metapopulation to population
ggplot(cv_diff, aes(x = combo_type, y = pop_cv_diff, group = Combo, fill = Combo)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
 # geom_errorbar(aes(ymin = phi_S_L2R_mean - phi_S_L2R_standard_error, ymax = phi_S_L2R_mean + phi_C_L2R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  labs(title = "Population to Metapopulation CV Dampening")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#########No summary of output ----- 
##Just extracting the raw values from partition results
partition_results_all <- do.call(rbind, lapply(partition_results, as.data.frame)) 
colnames(partition_results_all) <- gsub("\\.", "-", colnames(partition_results_all))  

partition_results_all <- partition_results_all %>%
  rownames_to_column() %>%
  separate(rowname, into = c("Iteration", "Metric"), sep = "\\.") 

partition_results_all <- partition_results_all %>%
  pivot_longer(cols = 3:12,
               names_to = "Combo",
               values_to = "Value") %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  mutate(combo_type = case_when(
    startsWith(Combo, "Glacier-fed_Glacier-fed_Glacier-fed") ~ "Homogenous",
    startsWith(Combo, "Glacier-fed_Glacier-fed_Rain-fed") ~ "2:1",
    startsWith(Combo, "Glacier-fed_Glacier-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Glacier-fed_Rain-fed_Rain-fed") ~ "2:1",
    startsWith(Combo, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combo, "Glacier-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combo, "Rain-fed_Rain-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combo, "Snow-fed_Snow-fed_Snow-fed") ~ "Homogenous",
  )) %>%
  mutate(cv_diff = CV_C_L - CV_C_R) %>%
  mutate(taxa = "inverts")

write.csv(partition_results_all, "data/variance_partitioning_results/invert_var_partition_results_bootstrapped.csv")

partition_results_all$combo_type <- ordered(partition_results_all$combo_type, 
                                            levels = c("Homogenous", "2:1", "Heterogenous"))
##Metacommunity stability
invert_meta_cv <- ggplot(partition_results_all, aes(x = combo_type, y = CV_C_R, group = combo_type, fill = combo_type)) +
  geom_boxplot()+
  #  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  ylab("Metacommunity Variability (CV-C,R)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
invert_meta_cv

##invert mean 
invert_meta_mean <- ggplot(partition_results_all, aes(x = combo_type, y = mean_metacom, group = combo_type, fill = combo_type)) +
  geom_boxplot()+
  #  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  ylab("Metacommunity Mean")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
invert_meta_mean

invert_meta_sd <- ggplot(partition_results_all, aes(x = combo_type, y = sd_metacom, group = combo_type, fill = combo_type)) +
  geom_boxplot()+
  #  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  ylab("Metacommunity SD")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
invert_meta_sd

##this looks fucked, something is not right.... 
##Local Community stability
invert_local_cv <- ggplot(partition_results_all, aes(x = combo_type, y = CV_C_L, group = combo_type, fill = combo_type)) +
  geom_boxplot()+
  #  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  ylab("Local Community Variability (CV-C,L) ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
invert_local_cv

invert_cv_diff <- ggplot(partition_results_all, aes(x = combo_type, y = cv_diff, group = combo_type, fill = combo_type)) +
  geom_boxplot()+
  #  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_standard_error, ymax = CV_C_R_mean + CV_C_R_standard_error), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  ylab("Local to Metacommunity Stabilization \n(Metacommunity CV-Local Community CV)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
invert_cv_diff
##
##Asynchrony -- community asynchrony local to regional -- driven by spatial community dissimilarity 
##Is the spatial synchrony of total community biomass among local patches
invert_comm_async <- ggplot(partition_results_all, aes(x = combo_type, y = phi_C_L2R, group = combo_type, fill = combo_type)) +
  geom_boxplot()+
  theme_classic() +
  ylab("Community Level Spatial Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
invert_comm_async

##library
library(ggpubr)

##Final plot 



invert_fig <- ggarrange(invert_meta_cv, invert_local_cv, invert_cv_diff, invert_comm_async, legend = "none", 
                          labels = c("a)", "b)", "c)", "d)"),
                          ncol = 2, nrow = 2, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
invert_fig

invert_results <- ggarrange()
sync_1_aov <- aov(phi_C_L2R ~ combo_type, data = partition_results_all)
summary(sync_1_aov)
TukeyHSD(sync_1_aov)

