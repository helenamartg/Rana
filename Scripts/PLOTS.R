#################################
# This script contains some plots
#################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

# Import phylogeny
library(phytools)
library(ape)
phylorana <- read.tree("DATA/phyloRana_74sp.tre")
phylorana$tip.label
phylorana$tip.label <- gsub("ana_", ". ", phylorana$tip.label)

# Import traits
dt74 <- read.csv("DATA/dt_74sp.csv")
dt74$sp <- gsub("ana_", ". ", dt74$sp)
rownames(dt74) <- dt74$sp

# Make dataset match with phylogeny
rownames(dt74) ==  phylorana$tip.label
dt74 <- dt74[phylorana$tip.label,]
svl <- setNames(dt74$SVL, rownames(dt74))
logsvl <- setNames(dt74$log_svl, rownames(dt74))

#plot svl across tree
contMap(phylorana, log(svl), fsize=.5, type='fan')

layout(1,1)
obj<-contMap(phylorana, svl, plot=FALSE)
plotTree.wBars(obj$tree, svl, fsize=0.4,scale=0.1,tip.labels=TRUE,
               type="fan",method="plotSimmap",colors=obj$cols)
add.color.bar(0.8,cols=obj$cols,title="trait value",obj$lims,digits=2,
              prompt=FALSE,x=0.9*par()$usr[1],y=0.9*par()$usr[4])

svl.object<-edge.widthMap(phylorana, svl)
plot(svl.object, max.width=0.8,color="dimgrey",
     legend="log(svl)" ,min.width=0.001)


# Another plot
plotTree.barplot(phylorana, svl)
layout(1,1)
plotTree.wBars(phylorana, svl, tip.labels = F, fsize = 0.5)

# TREE with colors regarding continents
plot(phylorana)
nodelabels(cex=0.6)

pal <-c("gold", "darkorchid")
dt74$sp == phylorana$tip.label


# Plot by real region and svl by dots
tree  <- paintSubTree(phylorana,node=112,state="2",stem=T)
cols<-c("darkorchid", "gold", "gold")
names(cols)<-c(1,2,3)
plotSimmap(tree,cols,lwd=2,pts=F, fsize=0.55, ftype="i")
tree<-paintSubTree(tree,node=79,state="3",stem=F)
layout(matrix(c(1,2),ncol=2, nrow=1))
layout.show(2)
par(mar=c(4,1,1,1))
plotSimmap(tree,cols,lwd=3,pts=F, fsize=0.55, ftype = "i")
plot(NULL, xlim=c(0,10), ylim=c(1,length(phylorana$tip.label)),
     ylab="", xlab="", bty = "n", xaxt="n", yaxt="n")
points(c(1:74)~svl, type = "p", pch = 21, cex = 0.8, bg=pal[as.factor(dt74$real_reg)])


# import subgenera
dt74[dt74$sp=="R. luteiventris",5] <- "WestNAM"
dt74[dt74$sp=="R. boylii", 5] <- "WestNAM"
dt74[dt74$sp=="R. sierrae", 5] <- "WestNAM"
dt74[dt74$sp=="R. muscosa", 5] <- "WestNAM"
dt74[dt74$sp=="R. aurora", 5] <- "WestNAM"
dt74[dt74$sp=="R. cascadae", 5] <- "WestNAM"

# plot phylogeny by state
states <- setNames(dt74$state, dt74$sp)
unique(states)

nam<-names(states)[states=="NAM"]
mex <- names(states)[states=="MEX"]
neo <- names(states)[states=="NEO"]
eu <- names(states)[states=="EU"]
as <- names(states)[states=="AS"]
westnam <-  names(states)[states=="WestNAM"]


library(tiff)
tiff("Fig1.tiff", compression = "lzw", width = 18, height = 18, units = "cm", res = 600)

tt<-paintBranches(phylorana,edge=sapply(as,match,phylorana$tip.label),
                  state="as", anc.state = "0")
tt<-paintBranches(tt,edge=sapply(nam,match,phylorana$tip.label),
                  state="nam")
tt<-paintBranches(tt,edge=sapply(westnam,match,phylorana$tip.label),
                  state="westnam")
tt<-paintBranches(tt,edge=sapply(mex,match,phylorana$tip.label),
                  state="mex")
tt <- paintBranches(tt,edge=sapply(neo,match,phylorana$tip.label),
                    state="neo")
tt <- paintBranches(tt,edge=sapply(eu,match,phylorana$tip.label),
                    state="eu")
cols<-setNames(c("grey68","darkorchid1","gold3","#735c08","gold4", "gold", "darkorchid4"),
               c("0","as","nam","westnam","mex", "neo", "eu"))
plot(tt,colors=cols,lwd=3,split.vertical=TRUE,ftype="i", fsize=0.6)

dev.off()

# Plot state With svl
library(tiff)
tiff("Fig1.tiff", compression = "lzw", width = 18, height = 18, units = "cm", res = 600)
layout(matrix(c(1,2),ncol=2, nrow=1))
par(mar=c(14, 4, 2, 2), xpd = T, family="sans", oma=c(2,2,0,0))
plot(tt,colors=cols,lwd=2.5,split.vertical=TRUE,ftype="i", fsize=0.6)

axis(1, at = seq(-3.088619, max(branching.times(phylorana)), length.out = 4),
     labels = round(seq(50, 0, length.out = 4), 1),
     cex.axis = 0.7, line = -0.7, padj=-1.5)

# axisPhylo(cex=0.1, line= -0.7)
mtext("MYA", side=1, line=0.5, at=23, cex=0.7)

barplot(svl, horiz = T, yaxt="n", xlim=c(0,500), space=0.5, xaxt="n")
axis(1, at = seq(-0, 200, length.out = 5), 
     labels = round(seq(0, 200, length.out = 5),2),
     cex.axis = 0.7, line = -0.7, padj=-1.5)
mtext("SVL (mm)", side=1, line=0.5, at=80, cex=0.7)
dev.off()


# Plot biomes across phylogeny
biomes <- setNames(dt74$biome, dt74$sp)
trop<-names(biomes)[biomes=="Tropical forest"]
temp <- names(biomes)[biomes=="Temperate Forest"]
gras <- names(biomes)[biomes=="Grasslands"]
boreal <- names(biomes)[biomes=="Boreal forest"]
des <- names(biomes)[biomes=="Desert"]

tt<-paintBranches(phylorana,edge=sapply(trop,match,phylorana$tip.label),
                  state="trop", anc.state = "0")
tt<-paintBranches(tt,edge=sapply(nam,match,phylorana$tip.label),
                  state="nam")
tt<-paintBranches(tt,edge=sapply(temp,match,phylorana$tip.label),
                  state="temp")
tt <- paintBranches(tt,edge=sapply(gras,match,phylorana$tip.label),
                    state="gras")
tt <- paintBranches(tt,edge=sapply(boreal,match,phylorana$tip.label),
                    state="boreal")
tt <- paintBranches(tt,edge=sapply(des,match,phylorana$tip.label),
                    state="des")
pal2 <- c("gray70","cadetblue3", "antiquewhite3", "chocolate", "darkolivegreen4", "chartreuse3")
cols<-setNames(pal2,
               c("0","boreal","des","gras", "temp", "trop"))
plot(tt,colors=cols,lwd=3,split.vertical=TRUE,ftype="i", fsize=0.6)



# plot svl across phylogeny
layout(1,1)
contMap(phylorana, svl, fsize=0.6)
plotBranchbyTrait(phylorana, svl)

# plot clutch size across phylogeny
clutch <- setNames(dt74$resid_clutch, dt74$sp)
plotBranchbyTrait(phylorana, clutch) # Not working because there are missing values

# plot ecotypes
eco <- setNames(dt74$ecotype, dt74$sp)
cols <- setNames(c("dodgerblue4", "cornflowerblue", "chocolate4"), unique(eco))
dotTree(phylorana, eco,colors = cols)
plotTree(phyloRana)

library(ggplot2)
library(ggtree)
state <- setNames(dt74$state, dt74$sp)
p <- ggtree(phylorana, aes(color=group), layout='circular') 

dt <- data.frame(id=phylorana$tip.label, value = dt74$log_svl)
p3 <- facet_plot(p, panel='bar', data=dt, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y), size=3, color='blue4') 
p3 + theme_tree2()


