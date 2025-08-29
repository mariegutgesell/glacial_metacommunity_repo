##Calculate temporal beta diversity for all invert bootstrapped communities 

##Load libraries
library(tidyverse)
library(reshape2)

##source beta diversity function
source("code/temporal_beta_diversity/lamy_spat_stability_function.R")

##import matrix in proper format from stacked array of bootstrapped invert communities 
source("code/temporal_beta_diversity/mixed_metacommunity_matrix.R")

source("code/sp_matrix_arrays/invert_sp_matrix_bootstrapped.R")

rm(list = ls()[!ls() %in% c("xtabs_list_collection")])


#convert spp matrix to community biomass in patch i at time t
sp_long<-matrix_format
s<- length(unique(sp_long$Site))

#number of sampled months (t)
t<-length(unique(sp_long$Time))

#community matrix (Y)

Y<-as.matrix(sp_long[ , -c(1, 2)])
str(Y)
#calculate metricsfrom spat_stab function####

res <- space_stab(Y, s, t)
mult.s <- melt(res)
str(mult.s)
attach(mult.s )

mult.s$variable<-as.factor(mult.s$variable)
mult.s$variable<-factor(mult.s$variable, levels=c("AlphaCV","GammaCV","PhiCV","AlphaHBD","GammaHBD","PhiHBD"))




##trying to get partitoning results for each of the 50 different arrays
tb_results <- list()

for(n in 1:50) {
  current_xtabs_list <- xtabs_list_collection[[paste0("Iteration_", n)]]
  
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
  mutate(taxa = "inverts")


tb_results_all_summary <- tb_results_all %>%
  group_by(combo_type) %>%
  summarise_at(vars(GammaHBD, PhiHBD), list(mean = mean, sd = sd))

tb_results_all_summary$combo_type <- ordered(tb_results_all_summary$combo_type, 
                                             levels = c("Homogenous", "2:1", "Heterogenous"))
##Plot out 
ggplot(tb_results_all_summary, aes(x = combo_type, y = mean, fill = combo_type)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean+sd))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(tb_results_all_summary, aes(x = combo_type, y = PhiHBD_mean, fill = combo_type)) +
  geom_col() +
  geom_errorbar(aes(ymin = PhiHBD_mean - PhiHBD_sd, ymax = PhiHBD_mean+PhiHBD_sd))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##put metacomm turnover and synchrony in 1 fig 
##do nmbds temporal plots for glacial gradient for supp 



library(vegan)
#plot out NMDS####
#### Multivariate dataset ####
# List of local communities
SiteL <- list(); for(i in 1:s) SiteL [[i]] <- Y[c(((i-1)*t+1):(i*t)),] 
# Metacommunity (sum across all local communities)
MetacomTime <- Reduce("+", SiteL) 
# hellinger transformation metacommunity
bio_h_meta <- decostand(MetacomTime, "hellinger")
# combined dataset of local communities and metacommunity
Y.all <- rbind(Y, MetacomTime)
bio.h.all <- decostand(Y.all, method="hell")

#### NMDS ####
nmds2k.all <- metaMDS(bio.h.all, distance="euclidean", trace=TRUE, trymax=1000, k=2)
nmds2k.all
stressplot(nmds2k.all, dist(bio.h.all)) # stress = 0.180
# extract scrs
sites_scrs <- as.data.frame(scores(nmds2k.all, display="sites"))
spps_scrs  <- as.data.frame(scores(nmds2k.all, display="species"))
spps_scrs$Species <- rownames(spps_scrs)
# compute axis ranges
xlim_scrs <- range(sites_scrs[,1], sites_scrs[,1])
ylim_scrs <- range(sites_scrs[,2], sites_scrs[,2])
#use plot info from NMDS script
# segments between local points

segments.loc <- data.frame()
for(j in 1:s) for(i in (t*(j-1)+1):(t*j-1)) segments.loc <- rbind(segments.loc, data.frame(x=sites_scrs[i,1], y=sites_scrs[i,2], xend=sites_scrs[i+1,1], yend=sites_scrs[i+1,2], col=col[j])) 
# start points
points.start <- data.frame()
for(j in 1:s) points.start <- rbind(points.start, data.frame(x=sites_scrs[t*(j-1)+1,1], y=sites_scrs[t*(j-1)+1,2], col=col[j]))
# add start point for the metacommunity
points.start <- rbind(points.start, data.frame(x=sites_scrs[t*3+1,1], y=sites_scrs[t*3+1,2], col="black"))#change '5' to whatever number of sites you have
points.start$Site <- as.factor(c(levels(matrix_format$Site), "Metacommunity"))
points.start$Site <- ordered(points.start$Site, levels=c(levels(matrix_format$Site), "Metacommunity"))
# end arrows
points.end <- data.frame()
for(j in 1:s) points.end <- rbind(points.end, data.frame(x=sites_scrs[t*j-1,1], y=sites_scrs[t*j-1,2], xend=sites_scrs[t*j,1], yend=sites_scrs[t*j,2], col=col[j])) 
# segments for the metacommunity
segments.reg <- data.frame()
for(i in (t*3+1):(t*5-1)) segments.reg <- rbind(segments.reg, data.frame(x=sites_scrs[i,1], y=sites_scrs[i,2], xend=sites_scrs[i+1,1], yend=sites_scrs[i+1,2], col="black"))#again change numbers to match number of sites +metacommunity
# actual plot
col<-c("#08519C", "darkseagreen", "aquamarine2","black")
pD <- 
  ggplot() +
  geom_point(data=points.start, aes(x=x, y=y, col=Site), size=5) +
  # add or remove species scores
  geom_text(data=spps_scrs, aes(x=NMDS1, y=NMDS2, label=Species), size=3, col="grey") +
  geom_segment(data=segments.loc, aes(x=x, y=y, xend=xend, yend=yend), col=segments.loc$col, linewidth=1.2, alpha=0.75) +
  geom_segment(data=points.end, aes(x=x, y=y, xend=xend, yend=yend), col=points.end$col, arrow=arrow(length=unit(0.55, "cm"), type="closed"), linewidth=0.75) +
  geom_segment(data=segments.reg, aes(x=x, y=y, xend=xend, yend=yend), col="black", linewidth=1.7) +
  geom_segment(data=points.end, aes(x=sites_scrs[t*4-1,1], y=sites_scrs[t*4-1,2], xend=sites_scrs[t*4,1], yend=sites_scrs[t*4,2]), col="black", #again change '6' to number of sites + metacommunity
               arrow=arrow(length=unit(0.5, "cm"), type="closed")) +
  xlab("Axis 1") + ylab("Axis 2") + 
  labs(title="(D) NMDS (stress = 0.13)") +
  coord_equal() +
  theme_bw()  + 
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=c(c(col, "black")))

pD

