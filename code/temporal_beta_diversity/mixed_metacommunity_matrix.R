##Formatting boostrapped arrays for all community combinations in matrix for temporal beta diversity analysis
library(tidyverse)
library(reshape2)


invert_sp_matrix <- readRDS("data/intermediate_data/invert_sp_matrix.rds")


# Convert xtabs to data frame
inv.matrix_data <- as.data.frame(as.table(invert_sp_matrix))

# Rename columns for clarity
colnames(inv.matrix_data) <- c("Species", "Time", "Site", "Biomass")

# Reshape data so that Species are in columns, with Time and Site as rows
matrix_format <- dcast(inv.matrix_data, Time + Site ~ Species, value.var = "Biomass")

# Sort matrix_format by Time and Site
matrix_format <- matrix_format[order(matrix_format$Site, matrix_format$Time), ]

# Convert to matrix, removing Time and Site columns if you need only species data
#sp_long<-as.matrix(matrix_format[ , -c(1, 2)])
#sp_long<-droplevels(matrix_format)
#attach(sp_long)
#Y<-sp_long
#Y<-as.matrix(Y)
#Y<-inv.matrix_list[[5]]
#s<-3
#t<-10
#str(Y)
#Y<-as.matrix(Y)

    
#sp_long
