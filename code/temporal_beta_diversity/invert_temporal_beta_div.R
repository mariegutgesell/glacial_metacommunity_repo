##Calculate temporal beta diversity for all invert bootstrapped communities 

##Load libraries
library(tidyverse)
library(reshape2)



##import matrix in proper format from stacked array of bootstrapped invert communities 
inverts <- readRDS("data/intermediate_data/invert_sp_matrix_stacked_array_bootstrapped.rds")

##source beta diversity function
source("code/temporal_beta_diversity/lamy_spat_stability_function.R")


##trying to get partitoning results for each of the 50 different arrays
tb_results <- list()

for(n in 1:50) {
  current_xtabs_list <- inverts[[paste0("Iteration_", n)]]
  
  iteration_results <- list()
  for(combo in names(current_xtabs_list)) {
    df <- current_xtabs_list[[combo]]
    df_long <- as.data.frame(as.table(df))
    
    # Rename columns for clarity
    colnames(df_long) <- c("Species", "Time", "Site", "Biomass")
    
    # Reshape data so that Species are in columns, with Time and Site as rows
    df_mat <- dcast(df_long, Time + Site ~ Species, value.var = "Biomass")
    
    # Sort matrix_format by Time and Site
    df_mat <- df_mat[order(df_mat$Site, df_mat$Time), ]
    s <- length(unique(df_mat$Site))
    t<-length(unique(df_mat$Time))
    Y <- as.matrix(df_mat[ , -c(1, 2)])
    tb_result <- space_stab(Y, s, t)
    iteration_results[[combo]] <- tb_result
  
  }
  
  iteration_results_df <- do.call(rbind, iteration_results)

  tb_results[[paste0("Iteration_", n)]] <- iteration_results_df
  
  
}


tb_results_all <- do.call(rbind, tb_results) 
colnames(tb_results_all) <- gsub("\\.", "-", colnames(tb_results_all))  

tb_results_all <- tb_results_all %>%
  rownames_to_column() %>%
  separate(rowname, into = c("Iteration", "Combo"), sep = "\\.") 

tb_results_all <- tb_results_all %>%
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
  mutate(HBD_diff = AlphaHBD - GammaHBD) %>%
  mutate(taxa = "inverts")


write.csv(tb_results_all, "data/temporal_betadiv_results/invert_temporal_betadiv_results_bootstrapped.csv")

tb_results_all_summary <- tb_results_all %>%
  group_by(taxa, combo_type) %>%
  summarise_at(vars(GammaCV, AlphaCV, PhiCV, GammaHBD, AlphaHBD, PhiHBD), list(mean = mean, sd = sd))

tb_results_all_summary$combo_type <- ordered(tb_results_all_summary$combo_type, 
                                             levels = c("Homogenous", "2:1", "Heterogenous"))

HBD_diff <- tb_results_all_summary %>%
  mutate(HBD_diff = AlphaHBD_mean - GammaHBD_mean)

HBD_diff$combo_type <- ordered(HBD_diff$combo_type, 
                                             levels = c("Homogenous", "2:1", "Heterogenous"))

##Plot out 
invert_mcv <- ggplot(tb_results_all_summary, aes(x = combo_type, y = GammaHBD_mean, fill = taxa, group = combo_type)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = GammaHBD_mean - GammaHBD_sd, ymax = GammaHBD_mean+GammaHBD_sd), position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue"))+
  scale_colour_manual(values = c("deepskyblue"))+
  ylab("Metacommunity Compositional\nVariability")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
invert_mcv

invert_mcs <- ggplot(tb_results_all_summary, aes(x = combo_type, y = PhiHBD_mean, fill = taxa, group = combo_type)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = PhiHBD_mean - PhiHBD_sd, ymax = PhiHBD_mean+PhiHBD_sd),  position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue"))+
  scale_colour_manual(values = c("deepskyblue"))+
  ylab("Spatial Compositional\nSynchrony")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
invert_mcs

invert_lcv <- ggplot(tb_results_all_summary, aes(x = combo_type, y = AlphaHBD_mean, fill = taxa, group = combo_type)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
  geom_errorbar(aes(ymin = AlphaHBD_mean - AlphaHBD_sd, ymax = AlphaHBD_mean+AlphaHBD_sd),  position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue"))+
  scale_colour_manual(values = c("deepskyblue"))+
  ylab("Local Community Compositional\nVariability")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
invert_lcv

##compositional variability dampening
invert_cvd <- ggplot(HBD_diff, aes(x = combo_type, y = AlphaHBD_mean, fill = taxa, group = combo_type)) +
  geom_col(aes(color = taxa), alpha = 0.5, linewidth =1) +
#  geom_errorbar(aes(ymin = AlphaHBD_mean - AlphaHBD_sd, ymax = AlphaHBD_mean+AlphaHBD_sd),  position = position_dodge(width = 0.9), width = 0.25)+
  theme_classic() +
  scale_fill_manual(values = c("deepskyblue"))+
  scale_colour_manual(values = c("deepskyblue"))+
  ylab("Compositional\nVariability Dampening")+
  theme(axis.title.y = element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.text = element_blank()) +
  facet_wrap(~taxa, nrow = 1)
invert_cvd

library(ggpubr)
invert_comp_fig <- ggarrange(invert_mcv, invert_lcv, invert_mcs, legend = "none", 
                             labels = c("a)", "b)", "c)"),
                             ncol = 1, nrow = 3, font.label = list(colour = "black", size = 14, family = "Times New Roman"))
invert_comp_fig


##summary statistics -- degree of decline in HBD 
0.2310398/0.1758152

##Glacial loss case study examples

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
  group_by(taxa, Combo, combo_type, streamID_2, Iteration, Month_Year, Family.x) %>%
  summarise_at(vars(Freq), list(iteration_total_biomass = sum)) %>%
  separate(streamID_2, into = c("stream_type", "rep"), sep = "_", remove =FALSE) %>%
  mutate(
    Family.x = as.character(Family.x),
    Month_Year = as.character(Month_Year),
    streamID_2 = as.character(streamID_2),
    iteration_total_biomass = as.numeric(iteration_total_biomass)
  ) %>%
  select(taxa, Combo, combo_type, streamID_2, Month_Year, rep, Family.x, iteration_total_biomass)
  


#inverts_df_total_biomass <- inverts_df %>%  
#  group_by(taxa, Combo, combo_type, Month_Year) %>%
#  summarise_at(vars(iteration_total_biomass), list(total_biomass = sum)) %>%
#  rename(Combination = "Combo")

#inverts_df_total_biomass$combo_type <- ordered(inverts_df_total_biomass$combo_type,
 #                                              levels = c("Heterogenous",  "2:1","Homogenous" )) 

inverts_df$combo_type <- ordered(inverts_df$combo_type,
                                 levels = c("Heterogenous",  "2:1","Homogenous" )) 



inverts_df$month <- months(as.Date(paste("01", inverts_df$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)
#inverts_df_total_biomass$month <- months(as.Date(paste("01", inverts_df_total_biomass$Month_Year, sep = "-"), format = "%d-%Y-%m"), abbr = TRUE)


#####################trying cleaned up code:  ###### --------
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(ggpubr)

## Initialize global axis limits
xlims_all <- c(Inf, -Inf)
ylims_all <- c(Inf, -Inf)

plot_nmds_trajectory <- function(df, combo_colors, label) {
  Y1 <- dcast(df, streamID_2 + Month_Year ~ Family.x, value.var = "iteration_total_biomass")
  s <- length(unique(Y1$streamID_2))
  t <- length(unique(Y1$Month_Year))
  Y1_num <- Y1[, 3:ncol(Y1)]
  
  # Remove empty species columns
  if (any(colSums(Y1_num, na.rm = TRUE) == 0)) {
    Y1_num <- Y1_num[, colSums(Y1_num, na.rm = TRUE) != 0]
  }
  
  if (nrow(Y1_num) == 0 || ncol(Y1_num) == 0) {
    stop("NMDS input data has no valid dimensions after removing 0-sum columns")
  }
  
  # Create local site list
  SiteL <- lapply(1:s, function(i) Y1_num[((i - 1) * t + 1):(i * t), , drop = FALSE])
  MetacomTime <- Reduce("+", SiteL)
  Y.all <- rbind(Y1_num, MetacomTime)
  bio.h.all <- decostand(Y.all, method = "hellinger")
  
  # Run NMDS
  nmds <- metaMDS(bio.h.all, distance = "euclidean", trymax = 1000, k = 2)
  
  if (any(is.na(scores(nmds, display = "sites")))) {
    stop("NMDS scores returned NA — likely failure to converge or invalid input matrix")
  }
  
  sites_scrs <- as.data.frame(scores(nmds, display = "sites"))
  spps_scrs <- as.data.frame(scores(nmds, display = "species"))
  spps_scrs$Species <- rownames(spps_scrs)
  
  # Update global axis limits
  assign("xlims_all", range(c(get("xlims_all", envir = .GlobalEnv), sites_scrs[, 1])), envir = .GlobalEnv)
  assign("ylims_all", range(c(get("ylims_all", envir = .GlobalEnv), sites_scrs[, 2])), envir = .GlobalEnv)
  
  # Create segments for local sites
  segments.loc <- data.frame()
  for (j in 1:s) {
    for (i in (t * (j - 1) + 1):(t * j - 1)) {
      segments.loc <- rbind(segments.loc, data.frame(
        x = sites_scrs[i, 1], y = sites_scrs[i, 2],
        xend = sites_scrs[i + 1, 1], yend = sites_scrs[i + 1, 2],
        col = combo_colors[j]
      ))
    }
  }
  
  # Start and end points
  points.start <- do.call(rbind, lapply(1:s, function(j)
    data.frame(x = sites_scrs[t * (j - 1) + 1, 1],
               y = sites_scrs[t * (j - 1) + 1, 2],
               col = combo_colors[j])))
  meta_start <- nrow(sites_scrs) - t + 1
  meta_end <- nrow(sites_scrs)
  
  points.start <- rbind(points.start,
                        data.frame(x = sites_scrs[meta_start, 1], y = sites_scrs[meta_start, 2], col = "black"))
  
  # Label site names
  Y1$streamID_2 <- as.factor(Y1$streamID_2)
  points.start$Site <- factor(c(levels(Y1$streamID_2), "Metacommunity"),
                              levels = c(levels(Y1$streamID_2), "Metacommunity"))
  
  points.end <- do.call(rbind, lapply(1:s, function(j)
    data.frame(x = sites_scrs[t * j - 1, 1],
               y = sites_scrs[t * j - 1, 2],
               xend = sites_scrs[t * j, 1],
               yend = sites_scrs[t * j, 2],
               col = combo_colors[j])))
  
  segments.reg <- data.frame()
  for (i in meta_start:(meta_end - 1)) {
    segments.reg <- rbind(segments.reg, data.frame(
      x = sites_scrs[i, 1], y = sites_scrs[i, 2],
      xend = sites_scrs[i + 1, 1], yend = sites_scrs[i + 1, 2],
      col = "black"
    ))
  }
  meta_arrow <- data.frame(
    x = sites_scrs[meta_end - 1, 1], y = sites_scrs[meta_end - 1, 2],
    xend = sites_scrs[meta_end, 1], yend = sites_scrs[meta_end, 2]
  )
  
  gg_theme <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.text = element_text(size = 8),
                    axis.title = element_text(size = 10),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                    legend.position = "bottom",
                    legend.title = element_blank())
  
  p <- ggplot() +
    geom_point(data = points.start, aes(x = x, y = y, col = Site), size = 5) +
    geom_segment(data = segments.loc, aes(x = x, y = y, xend = xend, yend = yend),
                 col = segments.loc$col, size = 0.75, alpha = 0.75) +
    geom_segment(data = points.end, aes(x = x, y = y, xend = xend, yend = yend),
                 col = points.end$col, arrow = arrow(length = unit(0.55, "cm"), type = "closed"), size = 0.75) +
    geom_segment(data = segments.reg, aes(x = x, y = y, xend = xend, yend = yend),
                 col = "black", size = 2) +
    geom_segment(data = meta_arrow, aes(x = x, y = y, xend = xend, yend = yend),
                 col = "black", arrow = arrow(length = unit(0.5, "cm"), type = "closed")) +
    xlab("Axis 1") + ylab("Axis 2") +
    coord_fixed() + theme_bw() + gg_theme +
    scale_color_manual(values = c(combo_colors, "black"))
  
  return(p)
}

### RUN FOR ALL THREE COMBINATIONS ###

# 1. Heterogeneous
df_het <- inverts_df %>% filter(Combo == "Glacier-fed_Rain-fed_Snow-fed")
pD_het <- plot_nmds_trajectory(df_het, c("deepskyblue4", "skyblue", "azure3"), "a)")
pD_het
# 2. 2:1 Snow-dominated
df_21 <- inverts_df %>% filter(Combo == "Rain-fed_Snow-fed_Snow-fed")
pD_21 <- plot_nmds_trajectory(df_21, c("skyblue", "azure3", "azure3"), "b)")
pD_21
# 3. Homogeneous
df_hom <- inverts_df %>% filter(Combo == "Rain-fed_Rain-fed_Rain-fed")
pD_hom <- plot_nmds_trajectory(df_hom, c("skyblue", "skyblue", "skyblue"), "c)")
pD_hom
# Apply unified axis limits to all plots
#pD_het <- pD_het + coord_cartesian(xlim = xlims_all, ylim = ylims_all)
#pD_21  <- pD_21  + coord_cartesian(xlim = xlims_all, ylim = ylims_all)
#pD_hom <- pD_hom + coord_cartesian(xlim = xlims_all, ylim = ylims_all)

# Combine into a panel
nmds_plot <- ggarrange(pD_het, pD_21, pD_hom, 
                       labels = c("a)", "b)", "c)"), 
                       ncol = 3, nrow = 1,
                       font.label = list(color = "black", size = 14, family = "Times New Roman"))

nmds_plot

