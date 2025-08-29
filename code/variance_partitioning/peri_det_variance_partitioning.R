##Partition Metacommunity Variation -- Periphyton and Detritus 


#Perohpyton/detritus:
source("code/sp_matrix_arrays/peri_det_sp_matrix.R")

rm(list = ls()[!ls() %in% c("peri_matrix", "det_matrix", "det_xtabs_list", "peri_xtabs_list")])



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
                        phi_S2C_R=phi_S2C_R)
  return(partition_3level)
}


##calculate CV and synchrony for periphyton
cv_peri <- var.partition(peri_matrix)
cv_peri <- as.data.frame(cv_peri)

##calculate CV and synchrony for detritus
cv_det <- var.partition(det_matrix)
cv_det <- as.data.frame(cv_det)



##Calculate for different site combinations 
cv_peri_2 <- lapply(peri_xtabs_list, var.partition)
peri_results_df <- do.call(rbind, lapply(names(cv_peri_2), function(name) {
  data <- as.data.frame(cv_peri_2[[name]])
  cbind(Combination = name, Metric = rownames(data), data)
})) %>%
  rename(Value = "cv_peri_2[[name]]") %>%
  spread(key = "Metric", value = "Value")


cv_det_2 <- lapply(det_xtabs_list, var.partition)
det_results_df <- do.call(rbind, lapply(names(cv_det_2), function(name) {
  data <- as.data.frame(cv_det_2[[name]])
  cbind(Combination = name, Metric = rownames(data), data)
})) %>%
  rename(Value = "cv_det_2[[name]]") %>%
  spread(key = "Metric", value = "Value")


##Spatial synchrony across organizational levels
ggplot(peri_results_df, aes(x = phi_S_L2R, y = phi_C_L2R)) +
  geom_point() +
  theme_classic() +
  xlab("Species-level spatial synchrony") +
  ylab("Community-level spatial synrchony")

##Species synchrony across scales
ggplot(peri_results_df, aes(x = phi_S2C_L, y = phi_S2C_R, color = Combination)) +
  geom_point() +
  theme_classic() +
  xlab("Local-scale species synchrony") +
  ylab("Regional-scale species synrchony")

##Plotting some results ----------
peri_results_df <- peri_results_df %>%
  mutate(combo_type = case_when(
    startsWith(Combination, "Glacier-fed_Glacier-fed_Glacier-fed") ~ "Homogenous",
    startsWith(Combination, "Glacier-fed_Glacier-fed_Rain-fed") ~ "2:1",
    startsWith(Combination, "Glacier-fed_Glacier-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Glacier-fed_Rain-fed_Rain-fed") ~ "2:1",
    startsWith(Combination, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combination, "Glacier-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combination, "Rain-fed_Rain-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Snow-fed_Snow-fed_Snow-fed") ~ "Homogenous",
  )) %>%
  mutate(taxa = "Periphyton")


det_results_df <- det_results_df %>%
  mutate(combo_type = case_when(
    startsWith(Combination, "Glacier-fed_Glacier-fed_Glacier-fed") ~ "Homogenous",
    startsWith(Combination, "Glacier-fed_Glacier-fed_Rain-fed") ~ "2:1",
    startsWith(Combination, "Glacier-fed_Glacier-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Glacier-fed_Rain-fed_Rain-fed") ~ "2:1",
    startsWith(Combination, "Glacier-fed_Rain-fed_Snow-fed") ~ "Heterogenous",
    startsWith(Combination, "Glacier-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Rain-fed_Rain-fed_Rain-fed") ~ "Homogenous",
    startsWith(Combination, "Rain-fed_Rain-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Rain-fed_Snow-fed_Snow-fed") ~ "2:1",
    startsWith(Combination, "Snow-fed_Snow-fed_Snow-fed") ~ "Homogenous",
  )) %>%
  mutate(taxa = "Detritus")

results_df <- rbind(peri_results_df, det_results_df) %>%
  mutate(cv_diff = CV_C_L - CV_C_R)

write.csv(results_df, "data/variance_partitioning_results/peri_det_var_partition_results.csv")

results_df$combo_type <- ordered(results_df$combo_type, 
        levels = c("Homogenous", "2:1", "Heterogenous"))

## PERIPHYTON - Metacommunity Variability
peri_meta_cv <- results_df %>%
  filter(taxa== "Periphyton") %>%
  ggplot(aes(x = combo_type, y = CV_C_R, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Metacommunity Variability (CV-C,R)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
peri_meta_cv

##Local community variability
peri_local_cv <- results_df %>%
  filter(taxa== "Periphyton") %>%
  ggplot(aes(x = combo_type, y = CV_C_L, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Local Community Variability (CV-C,L)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
peri_local_cv

peri_cv_diff <- results_df %>%
  filter(taxa== "Periphyton") %>%
  ggplot(aes(x = combo_type, y = cv_diff, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Local to Metacommunity Stabilization \n(Metacommunity CV-Local Community CV)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
peri_cv_diff

peri_comm_async <- results_df %>%
  filter(taxa== "Periphyton") %>%
  ggplot(aes(x = combo_type, y = phi_C_L2R, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Community Level Spatial Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
peri_comm_async

##final fig
peri_fig <- ggarrange(peri_meta_cv, peri_local_cv, peri_cv_diff, peri_comm_async, legend = "none", 
                        labels = c("a)", "b)", "c)", "d)"),
                        ncol = 2, nrow = 2, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
peri_fig


## DETRITUS- Metacommunity Variability
det_meta_cv <- results_df %>%
  filter(taxa== "Detritus") %>%
  ggplot(aes(x = combo_type, y = CV_C_R, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Metacommunity Variability (CV-C,R)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
det_meta_cv

##Local community variability
det_local_cv <- results_df %>%
  filter(taxa== "Detritus") %>%
  ggplot(aes(x = combo_type, y = CV_C_L, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Local Community Variability (CV-C,L)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
det_local_cv

det_cv_diff <- results_df %>%
  filter(taxa== "Detritus") %>%
  ggplot(aes(x = combo_type, y = cv_diff, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Local to Metacommunity Stabilization \n(Metacommunity CV-Local Community CV)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
det_cv_diff

det_comm_async <- results_df %>%
  filter(taxa== "Detritus") %>%
  ggplot(aes(x = combo_type, y = phi_C_L2R, group = combo_type, fill = combo_type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Community Level Spatial Synchrony")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
det_comm_async

##final fig
det_fig <- ggarrange(det_meta_cv, det_local_cv, det_cv_diff, det_comm_async, legend = "none", 
                      labels = c("a)", "b)", "c)", "d)"),
                      ncol = 2, nrow = 2, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
det_fig


