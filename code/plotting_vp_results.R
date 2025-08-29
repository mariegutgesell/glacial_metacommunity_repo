#### Plotting Variance Paritioning Results across all taxa

library(tidyverse)
library(ggplot2)
library(ggpubr)

##Read in variance partitoning results

inverts <- read.csv("data/variance_partitioning_results/invert_var_partition_results_bootstrapped.csv") %>%
  rename(Combination = "Combo") %>%
  select(taxa, Iteration, Combination, combo_type, CV_C_L, CV_C_R, phi_C_L2R, cv_diff, mean_metacom, sd_metacom)

pd <- read.csv("data/variance_partitioning_results/peri_var_partition_results_bootstrapped.csv") %>%
  rename(Combination = "Combo") %>%
  select(taxa, Iteration, Combination, combo_type, CV_C_L, CV_C_R, phi_C_L2R, cv_diff, mean_metacom, sd_metacom)

fish <- read.csv("data/variance_partitioning_results/fish_var_partition_results_bootstrapped.csv") %>%
  rename(Combination = "Combo") %>%
  select(taxa, Iteration, Combination, combo_type, CV_C_L, CV_C_R, phi_C_L2R, cv_diff, mean_metacom, sd_metacom)

df <- rbind(inverts, pd, fish) %>%
  mutate(CV_reduction_factor = CV_C_L/CV_C_R)

df$combo_type <- ordered(df$combo_type,
                         levels = c("Homogenous", "2:1", "Heterogenous"))
df$taxa <- ordered(df$taxa,
                   levels = c("Periphyton", "inverts", "fish"))

df_2 <- df 



##Plotting column chart of mean w/ standard deviation -- for all combinations ---------------------
standard_error <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
}

df_4 <- df_2 %>%
  # filter(Combination %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  select(taxa, Combination, combo_type, CV_C_R, CV_C_L, cv_diff, phi_C_L2R, mean_metacom, sd_metacom, CV_reduction_factor) %>%
  group_by(taxa, combo_type) %>%
  summarise_at(vars(CV_C_R, CV_C_L, cv_diff, phi_C_L2R, mean_metacom, sd_metacom, CV_reduction_factor), 
               list(mean = ~mean(.x, na.rm = TRUE), 
                    sd = ~sd(.x, na.rm = TRUE), 
                    se = ~standard_error(.x)))



df_4$combo_type <- ordered(df_4$combo_type,
                           levels = c("Homogenous", "2:1", "Heterogenous"))
df_4$taxa <- ordered(df_4$taxa,
                     levels = c("Periphyton", "inverts", "fish"))


meta_cv_2 <- ggplot(df_4, aes(x = combo_type, y = CV_C_R_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_sd, ymax = CV_C_R_mean + CV_C_R_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Metacommunity Biomass\nVariability (CV-C,R)")+
 # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
   facet_wrap(~taxa)
meta_cv_2


local_cv_2 <- ggplot(df_4, aes(x = combo_type, y = CV_C_L_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = CV_C_L_mean - CV_C_L_sd, ymax = CV_C_L_mean + CV_C_L_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Local Community Biomass\nVariability (CV-C,L)")+
 # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
    facet_wrap(~taxa)
local_cv_2



cv_diff_2<- ggplot(df_4, aes(x = combo_type, y = cv_diff_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = cv_diff_mean - cv_diff_sd, ymax = cv_diff_mean + cv_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  ylab("Local to Metacommunity Biomass Variability Dampening \n(Local Community CV-Metacommunity CV)")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
cv_diff_2


comm_async_2<- ggplot(df_4, aes(x = combo_type, y = phi_C_L2R_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = phi_C_L2R_mean - phi_C_L2R_sd, ymax = phi_C_L2R_mean + phi_C_L2R_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  ylab("Spatial Biomass\nSynchrony")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
comm_async_2

CV_rf<- ggplot(df_4, aes(x = combo_type, y = CV_reduction_factor_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = CV_reduction_factor_mean - CV_reduction_factor_sd, ymax = CV_reduction_factor_mean + CV_reduction_factor_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.8) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen", "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen", "deepskyblue", "darksalmon"))+
  ylab("CV Reduction Factor")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
CV_rf


meta_local_cv_2 <- ggarrange(meta_cv_2, local_cv_2, comm_async_2, legend = "none", 
                           labels = c("a)", "b)", "c)"),
                           ncol = 1, nrow = 3, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
meta_local_cv_2

final_fig <- ggarrange(meta_local_cv_2, cv_diff_2, legend = "none", ncol = 2, nrow = 1)
final_fig

###ANOVAs

meta_aov <- aov(CV_C_R ~ combo_type*taxa, df_2)
summary(meta_aov)
TukeyHSD(meta_aov)


##Plotting results for predicted climate homogenization gradient (rain, rain-rain-snow (or rain-snow-snow), and then rain-snow-ice) ------------

df_3 <- df_2 %>%
  filter(Combination %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  select(taxa, Combination, combo_type, CV_C_R, CV_C_L, cv_diff, phi_C_L2R) %>%
  group_by(taxa, Combination, combo_type) %>%
  summarise_at(vars(CV_C_R, CV_C_L, cv_diff, phi_C_L2R), 
                 list(mean = ~mean(.x, na.rm = TRUE), 
                      sd = ~sd(.x, na.rm = TRUE), 
                      se = ~standard_error(.x)))
    




df_3$combo_type <- ordered(df_3$combo_type,
                           levels = c("Heterogenous",  "2:1", "Homogenous"))
df_3$taxa <- ordered(df_3$taxa,
                     levels = c("Periphyton", "inverts", "fish"))


meta_cv_2 <- ggplot(df_3, aes(x = combo_type, y = CV_C_R_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = CV_C_R_mean - CV_C_R_sd, ymax = CV_C_R_mean + CV_C_R_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Metacommunity Biomass\nVariability (CV-C,R)")+
#  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
   facet_wrap(~taxa)
meta_cv_2


local_cv_2 <- ggplot(df_3, aes(x = combo_type, y = CV_C_L_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = CV_C_L_mean - CV_C_L_sd, ymax = CV_C_L_mean + CV_C_L_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Local Community Biomass\nVariability (CV-C,L)")+
 # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
  facet_wrap(~taxa)
local_cv_2


cv_diff_2<- ggplot(df_3, aes(x = combo_type, y = cv_diff_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = cv_diff_mean - cv_diff_sd, ymax = cv_diff_mean + cv_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  ylab("Local to Metacommunity Biomass Variability Dampening \n(Local Community CV-Metacmmunity CV)")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
 # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
  facet_wrap(~taxa, nrow = 1)
cv_diff_2


comm_async_2<- ggplot(df_3, aes(x = combo_type, y = phi_C_L2R_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = phi_C_L2R_mean - phi_C_L2R_sd, ymax = phi_C_L2R_mean + phi_C_L2R_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  ylab("Spatial Biomass\nSynchrony")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
 # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
   facet_wrap(~taxa, nrow = 1)
comm_async_2



meta_local_cv_2 <- ggarrange(meta_cv_2, local_cv_2, comm_async_2, legend = "none", 
                             labels = c("a)", "b)", "c)"),
                             ncol = 1, nrow = 3, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
meta_local_cv_2

final_fig <- ggarrange(meta_local_cv_2, cv_diff_2, legend = "none", ncol = 2, nrow = 1)
final_fig

cv_diff_peri <- df_3 %>%
  filter(taxa == "Periphyton") %>% 
  ggplot(aes(x = combo_type, y = cv_diff_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = cv_diff_mean - cv_diff_sd, ymax = cv_diff_mean + cv_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("darkgreen"))+
  scale_colour_manual(values = c("darkgreen"))+
  ylab("Local to Metacommunity Biomass\nVariability Dampening")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) #+
  # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
cv_diff_peri 


cv_diff_inverts <- df_3 %>%
  filter(taxa == "inverts") %>% 
  ggplot(aes(x = combo_type, y = cv_diff_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = cv_diff_mean - cv_diff_sd, ymax = cv_diff_mean + cv_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue"))+
  scale_colour_manual(values = c( "deepskyblue"))+
  ylab("Local to Metacommunity Biomass\nVariability Dampening")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) #+
# theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
cv_diff_inverts

cv_diff_fish <- df_3 %>%
  filter(taxa == "fish") %>% 
  ggplot(aes(x = combo_type, y = cv_diff_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = cv_diff_mean - cv_diff_sd, ymax = cv_diff_mean + cv_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("darksalmon"))+
  scale_colour_manual(values = c("darksalmon"))+
  ylab("Local to Metacommunity Biomass\nVariability Dampening")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) #+
# theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
cv_diff_fish

cv_diff_fig <- ggarrange(cv_diff_peri, cv_diff_inverts, cv_diff_fish, legend = "none", 
                         ncol = 1, nrow = 3, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
cv_diff_fig


#bring in time series fig
source("code/plot_combination_biomass.R")

fig_3 <- ggarrange(all_plot, cv_diff_fig, legend = "none", ncol = 2, nrow = 1, widths = c(2,1))
fig_3
##Plotting mean and SD of CV 
meta_mean <- ggplot(df_4, aes(x = combo_type, y = mean_metacom_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = mean_metacom_mean - mean_metacom_sd, ymax = mean_metacom_mean + mean_metacom_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Metacommunity Mean Biomass")+
   theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  #theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, scale = "free")
meta_mean

meta_sd <- ggplot(df_4, aes(x = combo_type, y = sd_metacom_mean, group = combo_type, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
   geom_errorbar(aes(ymin = sd_metacom_mean - sd_metacom_sd, ymax = sd_metacom_mean + sd_metacom_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("darkgreen",  "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Metacommunity Biomass SD")+
   theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  #theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, scale = "free")
meta_sd

sup_fig <- ggarrange(meta_mean, meta_sd, labels = c("a)", "b)"), legend = "none", ncol = 1, nrow = 2, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
sup_fig




