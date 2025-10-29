#### Plotting Temporal Beta Diversity Results across all taxa

library(tidyverse)
library(ggplot2)
library(ggpubr)

##Read in variance partitoning results

inverts <- read.csv("data/temporal_betadiv_results/invert_temporal_betadiv_results_bootstrapped.csv") %>%
  rename(Combination = "Combo") %>%
  select(taxa, Iteration, Combination, combo_type, GammaHBD, AlphaHBD, PhiHBD, HBD_diff)


fish <- read.csv("data/temporal_betadiv_results/fish_temporal_betadiv_results_bootstrapped.csv") %>%
  rename(Combination = "Combo") %>%
  select(taxa, Iteration, Combination, combo_type, GammaHBD, AlphaHBD, PhiHBD, HBD_diff)

df <- rbind(inverts,fish) %>%
  mutate(HBD_reduction_factor = AlphaHBD/GammaHBD) %>%
  mutate(combo_type_ab = case_when(
  startsWith(combo_type, "Het") ~ "Het",
  startsWith(combo_type, "2:1") ~ "2:1",
  startsWith(combo_type, "Hom") ~ "Hom",
))

df$combo_type_ab <- ordered(df$combo_type_ab,
                            levels = c("Hom", "2:1", "Het"))

df$taxa <- ordered(df$taxa,
                   levels = c("inverts", "fish"))

df_2 <- df 



##Plotting column chart of mean w/ standard deviation -- for all combinations ---------------------
standard_error <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
}

df_4 <- df_2 %>%
  # filter(Combination %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  select(taxa, Iteration, Combination, combo_type_ab, GammaHBD, AlphaHBD, PhiHBD, HBD_diff, HBD_reduction_factor) %>%
  group_by(taxa, combo_type_ab) %>%
  summarise_at(vars(GammaHBD, AlphaHBD, PhiHBD, HBD_diff, HBD_reduction_factor), 
               list(mean = ~mean(.x, na.rm = TRUE), 
                    sd = ~sd(.x, na.rm = TRUE), 
                    se = ~standard_error(.x)))



df_4$combo_type_ab <- ordered(df_4$combo_type_ab,
                           levels = c("Hom", "2:1", "Het"))
df_4$taxa <- ordered(df_4$taxa,
                     levels = c("inverts", "fish"))


meta_HBD_2 <- ggplot(df_4, aes(x = combo_type_ab, y = GammaHBD_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = GammaHBD_mean - GammaHBD_sd, ymax = GammaHBD_mean + GammaHBD_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Metacommunity Compositional\nVariability (GammaHBD)")+
  # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa)
meta_HBD_2


local_HBD_2 <- ggplot(df_4, aes(x = combo_type_ab, y = AlphaHBD_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = AlphaHBD_mean - AlphaHBD_sd, ymax = AlphaHBD_mean + AlphaHBD_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Community Compositional\nVariability (AlphaHBD)")+
  # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa)
local_HBD_2



HBD_diff_2<- ggplot(df_4, aes(x = combo_type_ab, y = HBD_diff_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = HBD_diff_mean - HBD_diff_sd, ymax = HBD_diff_mean + HBD_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("deepskyblue", "darksalmon"))+
  ylab("Compositional Variability Dampening \n(AlphaHBD-GammaHBD)")+
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
HBD_diff_2


comm_async_2<- ggplot(df_4, aes(x = combo_type_ab, y = PhiHBD_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = PhiHBD_mean - PhiHBD_sd, ymax = PhiHBD_mean + PhiHBD_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("deepskyblue", "darksalmon"))+
  ylab("Spatial Compositional\nSynchrony")+
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
comm_async_2


HBD_rf<- ggplot(df_4, aes(x = combo_type_ab, y = HBD_reduction_factor_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = HBD_reduction_factor_mean - HBD_reduction_factor_sd, ymax = HBD_reduction_factor_mean + HBD_reduction_factor_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.8) +
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c("deepskyblue", "darksalmon"))+
  ylab("HBD Reduction Factor")+
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
HBD_rf


meta_local_HBD_2 <- ggarrange(meta_HBD_2, local_HBD_2, comm_async_2, legend = "none", 
                             #labels = c("a)", "b)", "c)"),
                             ncol = 1, nrow = 3, font.label = list(colour = "black", size = 18, family = "Times New Roman"))
meta_local_HBD_2

final_fig <- ggarrange(meta_local_HBD_2, HBD_diff_2, legend = "none", ncol = 2, nrow = 1)
final_fig

###ANOVAs

meta_aov <- aov(CV_C_R ~ combo_type*taxa, df_2)
summary(meta_aov)
TukeyHSD(meta_aov)


##Plotting results for predicted climate homogenization gradient (rain, rain-rain-snow (or rain-snow-snow), and then rain-snow-ice) ------------

df_3 <- df_2 %>%
  filter(Combination %in% c("Glacier-fed_Rain-fed_Snow-fed", "Rain-fed_Rain-fed_Rain-fed", "Rain-fed_Snow-fed_Snow-fed")) %>%
  select(taxa, Combination, combo_type_ab, GammaHBD:HBD_reduction_factor) %>%
  group_by(taxa, Combination, combo_type_ab) %>%
  summarise_at(vars(GammaHBD, AlphaHBD, PhiHBD, HBD_diff, HBD_reduction_factor), 
               list(mean = ~mean(.x, na.rm = TRUE), 
                    sd = ~sd(.x, na.rm = TRUE), 
                    se = ~standard_error(.x)))




df_3$combo_type_ab <- ordered(df_3$combo_type_ab,
                           levels = c("Het",  "2:1", "Hom"))
df_3$taxa <- ordered(df_3$taxa,
                     levels = c("inverts", "fish"))


meta_tb_2 <- ggplot(df_3, aes(x = combo_type_ab, y = GammaHBD_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = GammaHBD_mean - GammaHBD_sd, ymax = GammaHBD_mean + GammaHBD_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c( "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c( "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Metacommunity Compositional\nVariability (GammaHBD)")+
  #  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
  facet_wrap(~taxa)
meta_tb_2


local_tb_2 <- ggplot(df_3, aes(x = combo_type_ab, y = AlphaHBD_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = AlphaHBD_mean - AlphaHBD_sd, ymax = AlphaHBD_mean + AlphaHBD_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c( "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c( "deepskyblue", "darksalmon"))+
  theme_classic() +
  ylab("Community Compositional\nVariability (AlphaHBD)")+
  # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
  facet_wrap(~taxa)
local_tb_2


tb_diff_2<- ggplot(df_3, aes(x = combo_type_ab, y =  HBD_diff_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin =  HBD_diff_mean -  HBD_diff_sd, ymax =  HBD_diff_mean +  HBD_diff_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c( "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c( "deepskyblue", "darksalmon"))+
  ylab("Compositional Variability Dampening \n(AlphaHBD-GammaHBD)")+
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
  facet_wrap(~taxa, nrow = 1)
tb_diff_2


comm_async_2<- ggplot(df_3, aes(x = combo_type_ab, y = PhiHBD_mean, group = combo_type_ab, fill = taxa)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = PhiHBD_mean - PhiHBD_sd, ymax = PhiHBD_mean + PhiHBD_sd, group = taxa), position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c( "deepskyblue", "darksalmon"))+
  scale_colour_manual(values = c( "deepskyblue", "darksalmon"))+
  ylab("Community Level \nSpatial Synchrony")+
  theme(axis.title.y = element_text(size = 18),axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  # theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  
  facet_wrap(~taxa, nrow = 1)
comm_async_2



meta_local_tb_2 <- ggarrange(meta_tb_2, local_tb_2, comm_async_2, legend = "none", 
                            
                             ncol = 1, nrow = 3, font.label = list(colour = "black", size = 18, family = "Times New Roman"))
meta_local_tb_2

final_fig <- ggarrange(meta_local_tb_2, tb_diff_2, legend = "none",  ncol = 2, nrow = 1, font.label = list(colour = "black", size = 18, family = "Times New Roman"))

final_fig




