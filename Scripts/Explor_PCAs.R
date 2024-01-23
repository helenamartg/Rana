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

# Add column with biomes
biomes <- read.csv("DATA/raw_data/biome_all.csv")
rownames(biomes) <- biomes$sp
biomes <- as.data.frame(biomes)
biomes$sp <- rownames(biomes)
biomes <- biomes[,5:6]

biom <- NULL
for (i in 1:nrow(biomes)){
  for (j in 1:nrow(morpho)){
    if (biomes[i,2]== morpho[j,2]){
      x <- biomes[i, c(1:2)]
      biom <- rbind(biom, x)
    }
  }
}
dim(biom)
dim(morpho)
morpho$sp==biom$sp
biom <- biom[morpho$sp,]
morpho <- cbind(morpho, biom[,1])


# Terrestrial vs aquatic
morpho[,20] <- NA
for (i in 1:nrow(morpho)){
  if (morpho[i,6] == "Terrestrial"){
    morpho[i,20] <- "Terrestrial"
  } else {
    morpho[i,20] <- "Aquatic"
  }
}



# PCA (without size)
morpho <- morpho[,c(2,4:6,14:17, 19,20)]
morpho.red <- na.omit(morpho) # 34 sp 
pca.glob <- prcomp(morpho.red[,5:8], scale. = T, center = T)
summary(pca.glob)
cor(morpho.red[,5:8], pca.glob$x)  

pal <- c("dodgerblue4", "cornflowerblue", "chocolate4")
plot(pca.glob$x[,1:2], bg=pal[as.factor(morpho.red$ecotype)], pch=21, cex=1.5)

pal0 <- c("dodgerblue3", "chocolate4")
plot(pca.glob$x[,1:2], bg=pal0[as.factor(morpho.red$V20)], pch=21, cex=1.5)

pal2 <- c("gold", "darkorchid")
plot(pca.glob$x[,1:2], col="white", bg=pal2[as.factor(morpho.red$phylo_reg)], 
     pch=21, cex=1.5, asp=1)

pal3 <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
plot(pca.glob$x[,1:2], col="white", bg=pal3[as.factor(morpho.red$state)], 
     pch=21, cex=1.5, asp=1)


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
p1 <- ggplot(data = morpho.red, aes(x = pc1, y = pc2, color=phylo_reg), scale= 0, asp=1) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal2) +
  theme_classic() 

# add polygons  
p1 + geom_mark_hull(concavity = 5,expand=0, radius=0, alpha=0.2, aes(fill=phylo_reg)) +
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


# PLOT phylomorphospace
#========================
# phylomorphospace: add phylogenetic relationships to a normal PCA

library(geomorph)
library(castor)

only.morph <- morpho.red[,5:8]
rownames(only.morph) <- morpho.red$sp

tr <- read.tree("DATA/phylorana_55sp.tre")
tr2 <- get_subtree_with_tips(tr, rownames(only.morph))$subtree

pca.morpho <- prcomp(only.morph, center. = T, scale. = T)
summary(pca.morpho)
cor(only.morph, pca.morpho$x)


layout(matrix(c(1,2), nrow=1, ncol=2))
pal2
pal2.alpha <- c(adjustcolor("gold", alpha.f = 0.6),
                adjustcolor("darkorchid", alpha.f = 0.6))


# By radiation
phylomorphospace(tr2, pca.morpho$x[,1:2], label="Off", node.size=c(0,0),
                 xlab = "", ylab = "", asp = 1, colors="grey40")
points(pca.morpho$x[,1], pca.morpho$x[,2], pch=21, 
       bg=pal2[as.factor(morpho.red$phylo_reg)], col="grey20",
       cex = 2)
mtext("PC1 - 79.16%", side=1, line=2.5, cex = 1.2, font = 2)
mtext("PC2 - 12.96%", side=2, line=2.5, cex = 1.2, font = 2)

plot(pca.morpho$x[,1], pca.morpho$x[,2], pch=21, 
     bg=pal2[as.factor(morpho.red$phylo_reg)], col="grey20",
     cex = 2, asp=1, ylab = "", xlab="")
phylomorphospace(tr2, pca.morpho$x[,1:2], label="Off", node.size=c(0,0),
                 xlab = "", ylab = "", colors="grey40", add=T)
mtext("PC1 - 79.16%", side=1, line=2.5, cex = 1.2, font = 2)
mtext("PC2 - 12.96%", side=2, line=2.5, cex = 1.2, font = 2)



# phyloPCA
pca.morpho
pca.w.phylo <- gm.prcomp(only.morph, phy = tr2, GLS = T)
summary(pca.w.phylo)
plot(pca.w.phylo, phylo = TRUE, asp=1, pch=21, col="grey20",
     bg=pal2[as.factor(morpho.red$phylo_reg)], 
     phylo.par = list(tip.labels = F, node.labels=F, node.cex=0),
     cex=1.5)


# ggplot
library(deeptime)
library(ggpmisc)

p1 <- ggplot(morpho.red, aes(x = pc1, y = pc2)) +
  geom_phylomorpho(tr2, color="grey40") +
  theme_classic() +
  xlab("PC1 (79.16%)") + 
  ylab("PC2 (12.96%)") 
  # scale_x_continuous(limits = symmetric_limits) +
  # scale_y_continuous(limits = symmetric_limits)

pca <- p1 + geom_point(size=3.5, aes(color=phylo_reg)) +
  scale_color_manual(values=pal2) +
  theme(legend.position="none")
  



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
p1 <- ggplot(data = morpho.ea, aes(x = pc1, y = pc2, color=V20), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal0[-c(3:5)]) +
  theme_classic() +
  ggtitle("Eurasia")

# add polygons  
ea <- p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=V20)) +
  scale_fill_manual(values=pal0[-c(3:5)]) +
  xlab("PC1 (81.54%)") + 
  ylab("PC2 (11.55%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))
# 
# # biomes
# pal6 <- c("cadetblue3","antiquewhite3", "chocolate", "darkolivegreen4", "chartreuse3")
# p1 <- ggplot(data = morpho.ea, aes(x = pc1, y = pc2, color=biom[,1]), scale= 0) +
#   geom_point(alpha=0.7, size=3) +
#   scale_color_manual(values=pal6) +
#   theme_classic() +
#   ggtitle("Eurasia")


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
p1 <- ggplot(data = morpho.am, aes(x = pc1, y = pc2, color=V19), scale= 0) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=pal0) +
  theme_classic() +
  ggtitle("America")

# add polygons  
am <- p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=V19)) +
  scale_fill_manual(values=pal0) +
  xlab("PC1 (74.59%)") + 
  ylab("PC2 (20.68%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=13),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=13))

plot_grid(ea, am)





# # Phylomorphospace by geographic region
#===========================================
# phylomorphospace(tr2, pca.morpho$x[,1:2], label="Off", node.size=c(0,0),
#                  xlab = "", ylab = "", asp = 1, colors=pal3)
# points(pca.morpho$x[,1], pca.morpho$x[,2], pch=21, col="grey20",
#        bg=pal3[as.factor(morpho.red$state)], 
#        cex = 2.3)
# mtext("PC1 - 79.16%", side=1, line=2.5, cex = 1.2, font = 2)
# mtext("PC2 - 12.96%", side=2, line=2.5, cex = 1.2, font = 2)



# PCA CLIMATIC
#==============
# clim <- read.csv("DATA/clim_74sp.csv")
# rownames(clim) <- clim$X
# clim <- clim[clim$X!="Rana_tlaloci",]
# clim <- clim[clim$X!="Rana_luteiventris",]
# clim <- clim[clim$X!="Rana_boylii",]
# clim <- clim[clim$X!="Rana_sierrae",]
# clim <- clim[clim$X!="Rana_muscosa",]
# clim <- clim[clim$X!="Rana_aurora",]
# clim <- clim[clim$X!="Rana_cascadae",]
# 
# clim <- clim[,-grep("sd", colnames(clim))]
# clim <- clim[,c(3:29)]
# 
# clim.var <- as.data.frame(scale(clim, center=T, scale=T))
# 
# pca.env <- prcomp(clim.var, scale. = T, center = T)
# summary(pca.env)
# cor(clim.var, pca.env$x)
# 
# # extract pc scores for first two component and add to dat dataframe
# clim.var$pc1 <- pca.env$x[,1]  # indexing the first column
# clim.var$pc2 <- pca.env$x[, 2]  # indexing the second column
# 
# pca.vars <- pca.env$rotation %>% data.frame
# pca.vars$vars <- rownames(pca.vars)
# pca.vars.m <- melt(pca.vars, id.vars = "vars")
# 
# # plot
# dt68 <- read.csv("DATA/dt_68sp.csv")
# rownames(dt68) <- dt68$sp
# rownames(clim.var) == dt68$sp 
# clim.var <- clim.var[dt68$sp,]
# 
# clim.var$state <- dt68$state 
# p1 <- ggplot(data = clim.var, aes(x = pc1, y = pc2, color=state), scale= 0) +
#   geom_point(alpha=0.7, size=3) +
#   scale_color_manual(values=pal3) +
#   theme_classic() 
# 
# # add polygons  
# p1 + geom_mark_hull(concavity = 5,expand=0,radius=0, alpha=0.2, aes(fill=state)) +
#   scale_fill_manual(values=pal3) +
#   xlab("PC1 (36.29%)") + 
#   ylab("PC2 (22.01%)") +
#   theme_classic() +
#   theme(axis.text.x = element_text(size=13),
#         axis.title = element_text(size=15),
#         axis.text.y = element_text(size=13))






