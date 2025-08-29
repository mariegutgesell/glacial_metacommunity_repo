
s<- length(unique(fishweight$Habitat))

#number of sampled months (t)
t<-length(unique(fishweight$month))

#community matrix (Y)
Y<-as.data.frame(community.s)#community.s is a matrix with time site and species 
which(apply(Y, 2, sum) == 0)

Y <- Y[,which(apply(Y,2, sum) != 0)]
Y<-as.matrix(Y)
#calculate metricsfrom spat_stab function####

res <- space_stab(Y, s, t)
mult.s <- melt(res)
str(mult.s)
attach(mult.s )

mult.s$variable<-as.factor(mult.s$variable)
mult.s$variable<-factor(mult.s$variable, levels=c("AlphaCV","GammaCV","PhiCV","AlphaHBD","GammaHBD","PhiHBD"))

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
points.start <- rbind(points.start, data.frame(x=sites_scrs[t*5+1,1], y=sites_scrs[t*5+1,2], col="black"))#change '5' to whatever number of sites you have
points.start$Site <- as.factor(c(levels(fishweight$Habitat), "Metacommunity"))
points.start$Site <- ordered(points.start$Site, levels=c(levels(fishweight$Habitat), "Metacommunity"))
# end arrows
points.end <- data.frame()
for(j in 1:s) points.end <- rbind(points.end, data.frame(x=sites_scrs[t*j-1,1], y=sites_scrs[t*j-1,2], xend=sites_scrs[t*j,1], yend=sites_scrs[t*j,2], col=col[j])) 
# segments for the metacommunity
segments.reg <- data.frame()
for(i in (t*5+1):(t*6-1)) segments.reg <- rbind(segments.reg, data.frame(x=sites_scrs[i,1], y=sites_scrs[i,2], xend=sites_scrs[i+1,1], yend=sites_scrs[i+1,2], col="black"))#again change numbers to match number of sites +metacommunity
# actual plot
pD <- 
  ggplot() +
  geom_point(data=points.start, aes(x=x, y=y, col=Site), size=5) +
  # add or remove species scores
  geom_text(data=spps_scrs, aes(x=NMDS1, y=NMDS2, label=Species), size=3, col="grey") +
  geom_segment(data=segments.loc, aes(x=x, y=y, xend=xend, yend=yend), col=segments.loc$col, linewidth=1.2, alpha=0.75) +
  geom_segment(data=points.end, aes(x=x, y=y, xend=xend, yend=yend), col=points.end$col, arrow=arrow(length=unit(0.55, "cm"), type="closed"), linewidth=0.75) +
  geom_segment(data=segments.reg, aes(x=x, y=y, xend=xend, yend=yend), col="black", linewidth=1.7) +
  geom_segment(data=points.end, aes(x=sites_scrs[t*6-1,1], y=sites_scrs[t*6-1,2], xend=sites_scrs[t*6,1], yend=sites_scrs[t*6,2]), col="black", #again change '6' to number of sites + metacommunity
               arrow=arrow(length=unit(0.5, "cm"), type="closed")) +
  xlab("Axis 1") + ylab("Axis 2") + 
  labs(title="(D) NMDS (stress = 0.16)") +
  coord_equal() +
  theme_bw() + gg_theme + 
  theme(legend.position="right", legend.title=element_blank()) +
  scale_color_manual(values=c(c(col, "black")))

pD

