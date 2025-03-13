##############################################
# Edit climatic variables: 
# remove 8 individuals without SVL max data.
###############################################

rm(list=ls())
library(phytools)
library(ape)
library(geomorph)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")


# Import climatic variables
#=============================
climEA <- read.csv("Data/clim_EA_31sp.csv")
climAM <- read.csv("Data/clim_AM_37sp.csv")

climEA <- climEA[,-c(1)]
rownames(climEA) <- climEA$X

climAM <- climAM[,-c(1)]
rownames(climAM) <- climAM$X


# Import trait dataset
dt <- read.csv("Data/dt_59sp.csv")
rownames(dt) <- dt$Species

dtEA <- dt[dt$phylo_subtree=="EA",]
dtAM <- dt[dt$phylo_subtree=="AM",]

climEA <- climEA[dtEA$Species,]
climAM <- climAM[dtAM$Species,]




# SAVE
# write.csv(climEA, "Data/climEA_26sp.csv")
# write.csv(climAM, "Data/climAM_33sp.csv")
