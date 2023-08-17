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


# PCA (without size)
morpho <- morpho[,c(2,4:6,14:17)]
morpho.red <- na.omit(morpho) # 34 sp 
pca.glob <- prcomp(morpho.red[,5:8], scale. = T, center = T)
summary(pca.glob)
cor(morpho.red[,5:8], pca.glob$x)  

pal <- c("dodgerblue4", "cornflowerblue", "chocolate4")
plot(pca.glob$x[,1:2], bg=pal[as.factor(morpho.red$ecotype)], pch=21, cex=1.5)

pal2 <- c("gold", "darkorchid")
plot(pca.glob$x[,1:2], col="white", bg=pal2[as.factor(morpho.red$phylo_reg)], pch=21, cex=1.5)

pal3 <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
plot(pca.glob$x[,1:2], col="white", bg=pal3[as.factor(morpho.red$state)], pch=21, cex=1.5)


##GGPLOT
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
p1 <- ggplot(data = morpho.red, aes(x = pc1, y = pc2, color=phylo_reg), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal2) +
  theme_classic() 

# add polygons  
p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=phylo_reg)) +
  scale_fill_manual(values=pal2) +
  xlab("PC1 (79.16%)") + 
  ylab("PC2 (12.06%)") +
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
  xlab("PC1 (79.16%)") + 
  ylab("PC2 (12.06%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))


# By ecotype
p3 <- ggplot(data = morpho.red, aes(x = pc1, y = pc2, color=ecotype), scale= 0) +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=pal) +
  theme_classic() 

p3 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=ecotype)) +
  scale_fill_manual(values=pal) +
  xlab("PC 1 (79.16%)") + 
  ylab("PC 2 (12.06%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))

# By subgen
subgen <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/subgen.csv", sep=";")
rownames(subgen) <- subgen$sp
subgen$sp == morpho.red$sp
subgen <- subgen[morpho.red$sp,]

morpho.red$subgen <- subgen$Subgenera
div <- morpho.red[morpho.red$sp!="Rana_weiningensis",]
div <- div[div$sp!="Rana_shuchinae",]
div <- div[div$subgen!="",]

pal5 <- c("#4F5D75", "#BFC0C0", "#9E0031", "#EF8354")

p3 <- ggplot(data = div, aes(x = pc1, y = pc2, color=subgen), scale= 0) +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=pal5) +
  theme_classic() 

p3 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=subgen)) +
  scale_fill_manual(values=pal5) +
  xlab("PC 1 (79.16%)") + 
  ylab("PC 2 (12.06%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))


# PCAs BY RADIATION
#===================
morpho
morpho.ea <- morpho[morpho$phylo_reg=="EA",]
morpho.am <- morpho[morpho$phylo_reg=="AM",]

### EURASIA
morpho.ea <- na.omit(morpho.ea)
pca.ea <- prcomp(morpho.ea[,5:8], scale. = T, center = T)
summary(pca.ea)
cor(morpho.ea[,5:8], pca.ea$x) 

# extract pc scores for first two component and add to dat dataframe
morpho.ea$pc1 <- pca.ea$x[,1]  # indexing the first column
morpho.ea$pc2 <- pca.ea$x[, 2]  # indexing the second column

pca.vars <- pca.ea$rotation %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# plot
p1 <- ggplot(data = morpho.ea, aes(x = pc1, y = pc2, color=state), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal3[-c(3:5)]) +
  theme_classic() 

# add polygons  
p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
  scale_fill_manual(values=pal3[-c(3:5)]) +
  xlab("PC1 (81.54%)") + 
  ylab("PC2 (11.55%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))


### AMERICA
morpho.am <- na.omit(morpho.am)
pca.am <- prcomp(morpho.am[,5:8], scale. = T, center = T)
summary(pca.am)
cor(morpho.am[,5:8], pca.am$x) 

# extract pc scores for first two component and add to dat dataframe
morpho.am$pc1 <- pca.am$x[,1]  # indexing the first column
morpho.am$pc2 <- pca.am$x[, 2]  # indexing the second column

pca.vars <- pca.am$rotation %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# plot
p1 <- ggplot(data = morpho.am, aes(x = pc1, y = pc2, color=state), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal3[-c(1,2)]) +
  theme_classic() 

# add polygons  
p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
  scale_fill_manual(values=pal3[-c(1,2)]) +
  xlab("PC1 (74.59%)") + 
  ylab("PC2 (20.68%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))


# PCA CLIMATIC
#==============
clim <- read.csv("DATA/clim_74sp.csv")
rownames(clim) <- clim$X
clim <- clim[clim$X!="Rana_tlaloci",]
clim <- clim[clim$X!="Rana_luteiventris",]
clim <- clim[clim$X!="Rana_boylii",]
clim <- clim[clim$X!="Rana_sierrae",]
clim <- clim[clim$X!="Rana_muscosa",]
clim <- clim[clim$X!="Rana_aurora",]
clim <- clim[clim$X!="Rana_cascadae",]

clim <- clim[,-grep("sd", colnames(clim))]
clim <- clim[,c(3:29)]

clim.var <- as.data.frame(scale(clim, center=T, scale=T))

pca.env <- prcomp(clim.var, scale. = T, center = T)
summary(pca.env)
cor(clim.var, pca.env$x)

# extract pc scores for first two component and add to dat dataframe
clim.var$pc1 <- pca.env$x[,1]  # indexing the first column
clim.var$pc2 <- pca.env$x[, 2]  # indexing the second column

pca.vars <- pca.env$rotation %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# plot
dt68 <- read.csv("DATA/dt_68sp.csv")
rownames(dt68) <- dt68$sp
rownames(clim.var) == dt68$sp 
clim.var <- clim.var[dt68$sp,]

clim.var$state <- dt68$state 
p1 <- ggplot(data = clim.var, aes(x = pc1, y = pc2, color=state), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal3) +
  theme_classic() 

# add polygons  
p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
  scale_fill_manual(values=pal3) +
  xlab("PC1 (36.29%)") + 
  ylab("PC2 (22.01%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))



library(geomorph)
# phylomorphospace: add phylogenetic relationships to a normal PCA
only.morph <- morpho.red[,5:8]
rownames(only.morph) <- morpho.red$sp

tr <- read.tree("DATA/phylorana_55sp.tre")
tr2 <- get_subtree_with_tips(tr, rownames(only.morph))$subtree

pca.morpho <- prcomp(only.morph, center. = T, scale. = T)
summary(pca.morpho)

# By radiation
phylomorphospace(tr2, pca.morpho$x[,1:2], label="Off", node.size=c(0,0),
                 xlab = "", ylab = "", asp = 1, colors=cols)
points(pca.morpho$x[,1], pca.morpho$x[,2], pch=21, 
       bg=pal2[as.factor(morpho.red$phylo_reg)], 
       cex = 1.6)
mtext("PC1 - 79.16%", side=1, line=2.5, cex = 1.2, font = 2)
mtext("PC2 - 12.96%", side=2, line=2.5, cex = 1.2, font = 2)

# By geographic region
phylomorphospace(tr2, pca.morpho$x[,1:2], label="Off", node.size=c(0,0),
                 xlab = "", ylab = "", asp = 1, colors=cols)
points(pca.morpho$x[,1], pca.morpho$x[,2], pch=21, 
       bg=pal3[as.factor(morpho.red$state)], 
       cex = 1.6)
mtext("PC1 - 79.16%", side=1, line=2.5, cex = 1.2, font = 2)
mtext("PC2 - 12.96%", side=2, line=2.5, cex = 1.2, font = 2)

