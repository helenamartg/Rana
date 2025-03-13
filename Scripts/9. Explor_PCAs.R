#################
# PCAs morpho
################

rm(list=ls())
library(phytools)
library(ape)
library(geomorph)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

dt<- read.csv("Data/dt_65sp.csv")
morpho <- read.csv("Data/morpho_phylores_45sp.csv")


# Extract 6 species from EA radiation that are distributed in aM
morpho <- morpho[morpho$Species!="Rana_tlaloci",]
morpho <- morpho[morpho$Species!="Rana_luteiventris",]
morpho <- morpho[morpho$Species!="Rana_boylii",]
morpho <- morpho[morpho$Species!="Rana_sierrae",]
morpho <- morpho[morpho$Species!="Rana_muscosa",]
morpho <- morpho[morpho$Species!="Rana_aurora",]
morpho <- morpho[morpho$Species!="Rana_cascadae",]


# PCA (without size)
morpho <- morpho[,c(2,4,13:16)]
morpho.red <- na.omit(morpho) # 26 sp 
pca.glob <- prcomp(morpho.red[,3:6], scale. = T, center = T)
summary(pca.glob)
cor(morpho.red[,3:6], pca.glob$x) 

pal0 <- c("darkorchid1", "darkorchid4")
pal1 <- c("gold4", "gold3", "gold")

pal2 <- c("gold", "darkorchid")
plot(pca.glob$x[,1:2], col="white", bg=pal2[as.factor(morpho.red$phylo_reg)], 
     pch=21, cex=1.5, asp=1)

pal3 <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
plot(pca.glob$x[,1:2], col="white", bg=pal3[as.factor(morpho.red$state)], 
     pch=21, cex=1.5, asp=1)


# GGPLOT
##########
library(ggplot2)
library(ggfortify)
library(cluster)

autoplot(pca.glob, scale=0, colour=pal3[as.factor(morpho.red$state)]) 

# extract pc scores for first two component and add to dat dataframe
morpho.red
morpho.red$pc1 <- pca.glob$x[,1]  # indexing the first column
morpho.red$pc2 <- pca.glob$x[, 2]  # indexing the second column

library(reshape2)
pca.vars <- pca.glob$rotation %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# By radiation
library(ggforce)
p1 <- ggplot(data = morpho.red, aes(x = pc1, y = pc2, color=phylo_reg), scale= 0, asp=1) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal2) +
  theme_classic() 

# add polygons  
p1 + geom_mark_hull(concavity = 5,expand=0, radius=0, alpha=0.2, aes(fill=phylo_reg)) +
  scale_fill_manual(values=pal2) +
  xlab("PC1 (82.37%)") + 
  ylab("PC2 (9.11%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))


# By state
p2 <- ggplot(data = morpho.red, aes(x = pc1, y = pc2, color=state), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal3) +
  theme_classic() 

p2 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
  scale_fill_manual(values=pal3) +
  xlab("PC1 (82.37%)") + 
  ylab("PC2 (9.11%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))


# PLOT phylomorphospace
#========================
# phylomorphospace: add phylogenetic relationships to a normal PCA

library(geomorph)
library(castor)

only.morph <- morpho.red[,3:6]
rownames(only.morph) <- morpho.red$Species

tr <- read.tree("Data/phylorana_65sp.tre")
tr2 <- get_subtree_with_tips(tr, rownames(only.morph))$subtree

pca.morpho <- prcomp(only.morph, center = T, scale. = T)
summary(pca.morpho)
cor(only.morph, pca.morpho$x)


layout(matrix(c(1,2), nrow=1, ncol=2))
pal2
pal2.alpha <- c(adjustcolor("gold", alpha.f = 0.6),
                adjustcolor("darkorchid", alpha.f = 0.6))

# 1. RAW PCA
#=============
# first phylomorphospace (the correct form)
par(mfrow=c(1,3), mar=c(4,4,15,2))
phylomorphospace(tr2, pca.morpho$x[,1:2], label="Off", node.size=c(0,0),
                 xlab = "", ylab = "", asp = 1, colors="grey60", ylim=c(-3,3))
points(pca.morpho$x[,1], pca.morpho$x[,2], pch=21, 
       bg=pal2[as.factor(morpho.red$phylo_reg)], col="grey20",
       cex = 1.5)
mtext("PC1 - 82.37%", side=1, line=2.5, cex = 1.2, font = 2)
mtext("PC2 - 9.11%", side=2, line=2.5, cex = 1.2, font = 2)
title(main="rawPCA (phylomorphospace)",
      font.main=3, cex=0.4)

# SAME
PCA.wphylo <- gm.prcomp(only.morph, phy = tr2)
summary(PCA.wphylo)
cor(only.morph, PCA.wphylo$x)
plot(PCA.wphylo,phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=1.5)
title(main="rawPCA + phylogenetic relationships",
      font.main=3, cex=0.4)

# 2. phyloPCA
#==============
# 2.1. PCA based on GLS-centering and projection (The one showed in Figure 1)
pca.w.phylo <- gm.prcomp(only.morph, phy = tr2, GLS = T)
summary(pca.w.phylo)
layout(1,1)
plot(pca.w.phylo, phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=2.5)
title(main="phyloPCA",
      font.main=3, cex=0.4)

eigenvalues <- pca.w.phylo$sdev^2
pca.w.phylo$rotation

# 2.2. PCA based on GLS-centering and transformed projection
pca.w.phylo.trans <- gm.prcomp(only.morph, phy = tr2, GLS = T, transform = T)
summary(pca.w.phylo.trans)
plot(pca.w.phylo.trans, phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=1.5)
title(main="phyloPCA transformed projection",
      font.main=3, cex=0.4)


# 3. PaCA
#==========
# 3.1. OLS method (rotation of PCA)
paca.ols <- gm.prcomp(only.morph, phy = tr2, align.to.phy = T)
plot(paca.ols, phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=1.5)
title(main="PACA - OLS method", font.main=3, cex=0.4)

# 3.2. GLS method (rotation of phylogenetic PCA)
paca.gls <- gm.prcomp(only.morph, phy = tr2, GLS = T, align.to.phy = T)
summary(paca.gls)
plot(paca.gls, phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=1.5)
title(main="PACA - GLS method", font.main=3, cex=0.4)

# 3.2. GLS method (rotation of phylogenetic PCA with transformation)
paca.gls.trans <- gm.prcomp(only.morph, phy = tr2, GLS = T, align.to.phy = T, transform = T)
plot(paca.gls.trans, phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=1.5)
title(main="PACA - GLS method and transformed projection", font.main=3, cex=0.4)



# ggplot
library(deeptime)
library(ggpmisc)

p1 <- ggplot(morpho.red, aes(x = pc1, y = pc2)) +
  geom_phylomorpho(tr2, color="grey40") +
  theme_classic() +
  xlab("PC1 (82.37%)") + 
  ylab("PC2 (9.11%)") 
  # scale_x_continuous(limits = symmetric_limits) +
  # scale_y_continuous(limits = symmetric_limits)

pca <- p1 + geom_point(size=3.5, aes(color=phylo_reg)) +
  scale_color_manual(values=pal2) +
  theme(legend.position="none")+
  ylim(-3,3)

par(mar=c(4,4,1,8))
pca

pca.morpho$sdev      # sdev of the PC (the square roots of the eigenvalues)
pca.morpho$rotation  # matrix whose columns contain the eigenvectors

# Calculate eigenvalues by squaring the standard deviations
eigenvalues <- pca.morpho$sdev^2
eigenvalues




# phyloPCAs BY RADIATION
#==========================
morpho
morpho.ea <- morpho[morpho$phylo_reg=="EA",]
morpho.am <- morpho[morpho$phylo_reg=="AM",]

### EURASIA
morpho.ea <- na.omit(morpho.ea)
rownames(morpho.ea) <- morpho.ea$Species
tr.ea <- (get_subtree_with_tips(tr2, morpho.ea$Species))$subtree

pca.EA.w.phylo <- gm.prcomp(morpho.ea[,3:6], phy = tr.ea, GLS = T)
summary(pca.EA.w.phylo)

pca.EA.w.phylo$sdev^2
pca.EA.w.phylo$rotation


# extract pc scores for first two component and add to dat dataframe
morpho.ea$pc1 <- pca.ea$x[,1]  # indexing the first column
morpho.ea$pc2 <- pca.ea$x[, 2]  # indexing the second column

pca.vars <- pca.ea$rotation %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# plot
p1 <- ggplot(data = morpho.ea, aes(x = pc1, y = pc2, color=state), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal0[-c(3:5)]) +
  theme_classic() +
  ggtitle("Eurasia")

# add polygons  
ea <- p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
  scale_fill_manual(values=pal0[-c(3:5)]) +
  xlab("PC1 (81.54%)") + 
  ylab("PC2 (11.55%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))

### AMERICA
morpho.am <- na.omit(morpho.am)
rownames(morpho.am) <- morpho.am$Species
tr.am <- (get_subtree_with_tips(tr2, morpho.am$Species))$subtree

pca.AM.w.phylo <- gm.prcomp(morpho.am[,3:6], phy = tr.am, GLS = T)
summary(pca.AM.w.phylo)

pca.AM.w.phylo$sdev^2
pca.AM.w.phylo$rotation

# extract pc scores for first two component and add to dat dataframe
morpho.am$pc1 <- pca.am$x[,1]  # indexing the first column
morpho.am$pc2 <- pca.am$x[, 2]  # indexing the second column

pca.vars <- pca.am$rotation %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# plot
p1 <- ggplot(data = morpho.am, aes(x = pc1, y = pc2, color=state), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal1) +
  theme_classic() +
  ggtitle("America")

# add polygons  
am <- p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
  scale_fill_manual(values=pal1) +
  xlab("PC1 (74.59%)") + 
  ylab("PC2 (20.68%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))

library(cowplot)
plot_grid(ea, am)






