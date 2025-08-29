# Getting Silhouttes of Animals

#Author(s): Marie
#Version: 2025-03-18

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)
library(rphylopic)
library(rsvg)
library(grid)

##### Code #####
#Create data frame for img coords
df <- data.frame(x = c(2, 4), y = c(10, 20))

##salmon
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Oncorhynchus nerka"), fill = "darksalmon",color="black", size = 10) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()

##dungeness crab
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Heptageniidae"), fill = "deepskyblue",color="black", size = 10) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()


##scallop
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Stylonematophyceae"), fill = "darkgreen",color="black", size = 10) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()
