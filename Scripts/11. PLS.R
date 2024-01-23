#############################################
# Morphology vs climatic variables (phyloPLS)
#############################################
# Only SVL 

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(phytools)
library(RRPP)
library(geiger)
library(pls)
library(geomorph)


# 1. Read phylogeny, morphometrics and environmental variables, and separate genera
#==================================================================================
phylorana_55 <- read.tree("DATA/phyloRana_55sp.tre")
phylorana_74 <- read.tree("DATA/phylorana_74sp.tre")
phylorana_68 <- read.tree("DATA/phylorana_68sp.tre")

clim <- read.csv("DATA/clim_74sp.csv")
rownames(clim) <- clim$X

traits <- read.csv("DATA/dt_morpho_55sp.csv")
rownames(traits) <- traits$sp
morpho <- traits[,c(13:17)]

dt <- read.csv("DATA/dt_74sp.csv")
rownames(dt) <- dt$sp


# 68 species
dt68 <- read.csv("DATA/dt_68sp.csv")
rownames(dt68) <-  dt68$sp

dt68 <- dt68[phylorana_68$tip.label,]
dt68$sp == phylorana_68$tip.label

clim68 <- clim[dt68$sp,]
clim68$X == dt68$sp


# 2. Check that phenotypic data and tree tips are in the same order
#==================================================================
rownames(morpho) ==  phylorana_55$tip.label
morpho <- morpho[phylorana_55$tip.label,]

rownames(dt) == phylorana_74$tip.label
dt <- dt[phylorana_74$tip.label,]

rownames(clim) == phylorana_55$tip.label
clim_55 <- clim[phylorana_55$tip.label,]
clim_55 <- clim_55[,-grep("sd", colnames(clim))]
clim_55 <- clim_55[,c(3:29)]

clim_74 <- clim[phylorana_74$tip.label,]
clim_74 <- clim_74[,-grep("sd", colnames(clim))]
clim_74 <- clim_74[,c(3:29)]

clim68 <- clim68[,-grep("sd", colnames(clim68))]
clim68 <- clim68[,c(3:29)]


## Scale climatic dataset for ALL
clim_55.var <- as.data.frame(scale(clim_55, center=T, scale = T))
clim_55.var <- as.matrix(clim_55.var)

clim_74.var <- as.data.frame(scale(clim_74, center=T, scale=T))
clim_74.var <- as.matrix(clim_74.var)

clim68.var <- as.data.frame(scale(clim68, center=T, scale=T))
clim68.var <- as.matrix(clim68.var)


# 3. phyloPLS
#==============
# size
svl <- as.matrix(setNames(dt68$log_svl, rownames(dt68)))

sz.phyloPLS <- phylo.integration(svl, clim68.var, phylorana_68)
sz.phyloPLS
summary(sz.phyloPLS)
plot(sz.phyloPLS, pch=21, bg="black", lwd=2)
plot(sz.phyloPLS$YScores ~ sz.phyloPLS$XScores, pch=21, bg="black")
abline(lm(sz.phyloPLS$YScores ~ sz.phyloPLS$XScores))


# shape (morpho)
library(castor)
shape <- as.matrix(morpho[,-1])
morpho2 <- na.omit(shape) 
tree2 <- get_subtree_with_tips(phylorana_55, rownames(morpho2))$subtree
clim_55.var <- clim_55.var[rownames(morpho2),]
rownames(clim_55.var) == rownames(morpho2)

sh.phyloPLS <- phylo.integration(morpho2, clim_55.var, tree2)
sh.phyloPLS
plot(sh.phyloPLS, pch=21, bg="black", lwd=2)

par(mar=c(6,4,3,3))
barplot(sz.phyloPLS$right.pls.vectors[,1], las = 2, font=2,  horiz = T)

rownames(sh.phyloPLS$left.pls.vectors) <- gsub("res_", "", rownames(sh.phyloPLS$left.pls.vectors))
barplot(sh.phyloPLS$left.pls.vectors[,1], las=2, font=2, ylim = c(-1,1))
barplot(sh.phyloPLS$right.pls.vectors[,1],las = 2, font=2, ylim=c(-0.6,0.6))



# 4. By radiation
#==================
# remove 6 species
traits
traits <- traits[traits$sp!="Rana_boylii",]
traits <- traits[traits$sp!="Rana_sierrae",]
traits <- traits[traits$sp!="Rana_muscosa",]
traits <- traits[traits$sp!="Rana_tlaloci",] # No climatic data

dt
dt <-  dt[dt$sp!="Rana_boylii",]
dt <- dt[dt$sp!="Rana_sierrae",]
dt <- dt[dt$sp!="Rana_muscosa",]
dt <- dt[dt$sp!="Rana_tlaloci",] # No climatic data


# Shape
am <- traits[traits$phylo_reg=="AM",]
morpho.am <- am[,c(14:17)]
morpho.am <- na.omit(morpho.am)

rownames(clim_55.var) <- as.factor(rownames(clim_55.var))
clim.am.sh <- clim_55.var[rownames(morpho.am),]
rownames(clim.am.sh) == rownames(morpho.am)

phylo.am.sh <- get_subtree_with_tips(phylorana_55, rownames(morpho.am))$subtree


ea <- traits[traits$phylo_reg=="EA",]
morpho.ea <- ea[,c(14:17)]
morpho.ea <- na.omit(morpho.ea)

clim.ea.sh <- clim_55.var[rownames(morpho.ea),]
rownames(clim.ea.sh) == rownames(morpho.ea)

phylo.ea.sh <- get_subtree_with_tips(phylorana_55, rownames(morpho.ea))$subtree


# size
svl.am <- dt[dt$phylo_reg=="AM",]
svl.am <- setNames(svl.am$log_svl, svl.am$sp)

clim.am.sz <- clim_74.var[names(svl.am),]
rownames(clim.am.sz) == names(svl.am)

phylo.am.sz <- get_subtree_with_tips(phylorana_74, names(svl.am))$subtree


svl.ea <- dt[dt$phylo_reg=="EA",]
svl.ea <- setNames(svl.ea$log_svl, svl.ea$sp)

clim.ea.sz <- clim_74.var[names(svl.ea),]
rownames(clim.ea.sz) == names(svl.ea)

phylo.ea.sz <- get_subtree_with_tips(phylorana_74, names(svl.ea))$subtree



# PLS EURASIA
#==============
sz.phyloPLS.ea <- phylo.integration(svl.ea, clim.ea.sz, phylo.ea.sz)
sz.phyloPLS.ea

sh.phyloPLS.ea <- phylo.integration(morpho.ea, clim.ea.sh, phylo.ea.sh)
sh.phyloPLS.ea

par(mfrow = c(1,2))
plot(sz.phyloPLS.ea, pch=21, col="black", bg="darkorchid", lwd=2, cex=1.8)
plot(sh.phyloPLS.ea, pch=21, col="black", bg="darkorchid", lwd=2, cex=1.8)


# AMERICA
sz.phyloPLS.am <- phylo.integration(svl.am, clim.am.sz, phylo.am.sz)
sz.phyloPLS.am

sh.phyloPLS.am <- phylo.integration(morpho.am, clim.am.sh, phylo.am.sh)
sh.phyloPLS.am


par(mfrow = c(1,2))
plot(sz.phyloPLS.am, pch=21, col="black", bg="gold", lwd=2, cex=1.8)
plot(sh.phyloPLS.am, pch=21, col="black", bg="gold", lwd=2, cex=1.8)


# PLOT ALL
par(mfrow = c(2,2))
plot(sz.phyloPLS.ea, pch=21, col="black", bg="darkorchid", lwd=2, cex=1.8)
plot(sh.phyloPLS.ea, pch=21, col="black", bg="darkorchid", lwd=2, cex=1.8)
plot(sz.phyloPLS.am, pch=21, col="black", bg="gold", lwd=2, cex=1.8)
plot(sh.phyloPLS.am, pch=21, col="black", bg="gold", lwd=2, cex=1.8)


# SAVE #
#### SAVE ####
sizes <- cbind(sz.phyloPLS$right.pls.vectors, sz.phyloPLS.ea$right.pls.vectors, 
               sz.phyloPLS.am$right.pls.vectors)
shapes.clim <- cbind(sh.phyloPLS$right.pls.vectors[,1], sh.phyloPLS.ea$right.pls.vectors[,1],
                     sh.phyloPLS.am$right.pls.vectors[,1])
shapes.morpho <- cbind(sh.phyloPLS$left.pls.vectors[,1], sh.phyloPLS.ea$left.pls.vectors[,1],
                       sh.phyloPLS.am$left.pls.vectors[,1])
shapes <- rbind(shapes.clim, shapes.morpho)

colnames(sizes) <- c("Global", "Eurasia", "America")
colnames(shapes) <- colnames(sizes)

# write.table(sizes, "Results/PLS_sizes_all.txt", sep=";")
# write.table(shapes, "Results/PLS_shapes_all.txt", sep=";")


