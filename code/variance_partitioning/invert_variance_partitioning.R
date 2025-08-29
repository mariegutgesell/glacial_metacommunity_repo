##Partition Metacommunity Variation -- Benthic inverts, non-bootstrapped samples 


##Inverts:
source("code/sp_matrix_arrays/invert_sp_matrix.R")


rm(list = ls()[!ls() %in% c("sp_matrix_2", "sp_matrix_all", "xtabs_list")])

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


##calculate CV and synchrony for all sites
cv_all <- var.partition(sp_matrix_all)

cv_all <- as.data.frame(cv_all) %>%
  mutate(value_squared = cv_all^2)


##calculate CV and synchrony for only 3 sites (excluding mixed)
cv_2 <- var.partition(sp_matrix_2)

cv_2 <- as.data.frame(cv_2)

##calculate CV and synchrony for each of the site combinations
cv_3 <- lapply(xtabs_list, var.partition)
results_df <- do.call(rbind, lapply(names(cv_3), function(name) {
  data <- as.data.frame(cv_3[[name]])
  cbind(Combination = name, Metric = rownames(data), data)
})) %>%
  rename(Value = "cv_3[[name]]") %>%
  spread(key = "Metric", value = "Value")



##Spatial synchrony across organizational levels
ggplot(results_df, aes(x = phi_S_L2R, y = phi_C_L2R)) +
  geom_point() +
  theme_classic() +
  xlab("Species-level spatial synchrony") +
  ylab("Community-level spatial synrchony")

##Species synchrony across scales
ggplot(results_df, aes(x = phi_S2C_L, y = phi_S2C_R, color = Combination)) +
  geom_point() +
  theme_classic() +
  xlab("Local-scale species synchrony") +
  ylab("Regional-scale species synrchony")

