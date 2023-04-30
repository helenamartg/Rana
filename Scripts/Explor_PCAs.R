#################
# PCAs morpho
################

rm(list=ls())
library(phytools)
library(ape)
library(geomorph)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

dt<- read.csv("DATA/dt_74sp.csv")
morpho <- read.csv("DATA/dt_morpho_55sp.csv")


# Extract 6 species from EA radiation that are distributed in aM
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]

# Reduce the dataset to species that have values for all morphological variables (same sample size)
morpho.red <- na.omit(morpho[,c(2:5,7:12)])
table(morpho.red$phylo_reg)

# Calculate Residuals
mat.red <- NULL
for (i in 7:ncol(morpho.red)){
  res <- resid(lm(log(morpho.red[,i]) ~ log(morpho.red$SVL)))
  mat.red <- cbind(mat.red, res)
}

colnames(mat.red) <- c("res_HL", "res_HW", "res_FL", "res_TL")
rownames(mat.red) <- morpho.red$sp

# Merge datasets
mat.red <- cbind(morpho.red, mat.red) # 34 sp

# PCA (without size)
pca.glob <- prcomp(mat.red[,11:14], scale. = T, center = T)
cor(mat.red[,11:14], pca.glob$x)  

pal <- c("dodgerblue4", "cornflowerblue", "chocolate4")
plot(pca.glob$x[,1:2], bg=pal[as.factor(mat.red$ecotype)], pch=21, cex=1.5)

pal2 <- c("gold", "darkorchid")
plot(pca.glob$x[,1:2], bg=pal2[as.factor(mat.red$phylo_reg)], pch=21, cex=1.5)

pal3 <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
plot(pca.glob$x[,1:2], bg=pal3[as.factor(mat.red$state)], pch=21, cex=1.5)

# By radiation
mat.red.EA <- mat.red[mat.red$phylo_reg=="EA",]
mat.red.AM <- mat.red[mat.red$phylo_reg=="AM",]

pca.ea <- prcomp(mat.red.EA[,6:10], scale. = T, center = T)
cor(mat.red.EA[,11:14], pca.ea$x)  

pal <- c("dodgerblue4", "cornflowerblue", "chocolate4")
plot(pca.ea$x[,1:2], bg=pal[as.factor(mat.red.EA$ecotype)], pch=21, cex=1.5)



# PCA (without size)
pca.glob <- prcomp(mat.red[,6:10], scale. = T, center = T)
cor(mat.red[,11:14], pca.glob$x)  

pal <- c("dodgerblue4", "cornflowerblue", "chocolate4")
plot(pca.glob$x[,1:2], bg=pal[as.factor(mat.red$ecotype)], pch=21, cex=1.5)

pal2 <- c("gold", "darkorchid")
plot(pca.glob$x[,1:2], bg=pal2[as.factor(mat.red$phylo_reg)], pch=21, cex=1.5)

pal3 <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
plot(pca.glob$x[,1:2], bg=pal3[as.factor(mat.red$state)], pch=21, cex=1.5)

