############################################
#   Testing diversification patterns  (HMG)
############################################

# This scripts contains:
#  1. Lineage-through-time plots (LTT)
#  2. Disparity-through-time plots (DTT)
#  3. Fitting different models of diversification

rm(list=ls())

library(phytools)
library(ape)
library(geomorph)
library(phyloch)
library(strap)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

# 1. Import phylogeny
#=====================
tr74 <- read.tree("DATA/phyloRana_74sp.tre")

# Original phylogeny
tr <- read.nexus("DATA/raw_data/Yuan_2016_phylo_v2.tree")
is.ultrametric(tr)
par(mar=c(1,1,1,1))
plot(tr, cex=0.6)
nodelabels()

# Remove outgroups
tr <- extract.clade(tr, 97)
plot(tr) #New tree has 90 tips

# Remove R_kunyuensis (= R. coreana)
tr$tip.label
tr <- drop.tip(tr, 2) #New tree with 89 tips
plot(tr, cex=.5)
is.ultrametric(tr)

# Separate in two radiations
plot(tr, cex=.5)
nodelabels(cex=0.5)

# tr <- rotate(tr, 91)
# tr <- root(tr, "Rana_weiningensis", resolve.root=T)
# plot(tr, cex=0.5)

phyloAM <- extract.clade(tr, 92)
plot(phyloAM, cex=0.7)  # 47 sp (8 undescribed species)

tr$tip.label
phyloEA <- drop.clade(tr, c(93:137))
plot(phyloEA, cex=0.7)
phyloEA$tip.label
phyloEA <- drop.tip(phyloEA, "NA")
plot(phyloEA)

# Remove 6 species
phyloEA$tip.label
nodelabels()  # Clade 46
phyloEA <- drop.clade(phyloEA, c(47:50))
phyloEA <- drop.tip(phyloEA, "NA")
plot(phyloEA) # 36 tips

is.ultrametric(phyloAM)
is.ultrametric(phyloEA)

# 2. Import traits
#===================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
dt74 <- dt74[,-1]
rownames(dt74) <- dt74$sp

# Transform SVL --> log(SVL)
dt74$logSVL <- log(dt74$SVL)

svl <- setNames(dt74$logSVL, dt74$sp)


# 3. LINEAGE-THROUGH-TIME
#=========================
par(mar=c(4,4,4,4))
ltt(phyloEA, main="European radiation")
ltt(phyloAM, main="American radiation")

par(mar=c(4,4,4,4))
plot(ltt(phyloEA), show.tree = T, main="Eurasian radiation", ylim=c(0,4))
plot(ltt(phyloAM), show.tree = T, main="American radiation", ylim=c(0,4))

# gamma statistic though time
gtt_EA <- gtt(phyloEA) 
plot(gtt_EA)
gtt_AM <- gtt(phyloAM)
plot(gtt_AM)


#4. The gamma statistic
#=======================
ltt_EA <- ltt(phyloEA)
ltt_EA$gamma
EA.mccr <- mccr(ltt_EA) #Takes into account incomplete taxon sampling in computing pval for gamma statistic
EA.mccr
plot(EA.mccr)

ltt_AM <- ltt(phyloAM)
ltt_AM$gamma
AM.mccr <- mccr(ltt_AM)
AM.mccr
plot(AM.mccr)

# PLOT BOTH LTT in one single plot
plot(ltt_AM$times, ltt_AM$ltt, type="p", pch=16, bg="black", lwd=2, 
     main="Lineage-through-time with all sp", xlab="time", ylab="nº of lineages)", xlim=c(0,50))
lines(ltt_EA$times, ltt_EA$ltt, type="p", pch=21, bg="darkgrey", lwd=2)
legend("bottomright", legend = c("American radiation", "European radiation"), 
       pt.bg = c("black", "darkgrey"), bty="n", pch=21)


# AXES Changed
plot(ltt_EA$times, ltt_EA$ltt, axes=F, type = "l", lwd=2, col="darkgrey", 
     ylim=c(0,50), xlab="Million years", ylab="nº of lineages")
axis(2, las=2,cex.axis=0.8)
labs<-axTicks(1)
h<-max(nodeHeights(tr))
at<-h-labs
axis(1,at=at,labels=labs,cex.axis=0.8)

# Transform times of American radiation to plot together
times2 <- NULL
for (i in 1:length(ltt_AM$times)){
  times2[i] <- ltt_AM$times[i] + 9.310806
}
lines(times2, ltt_AM$ltt, lwd=2, col="black")
box()
legend("topleft", legend = c("American radiation", "European radiation"), 
       pt.bg = c("black", "darkgrey"), bty="n", pch=21)



# 5. DISPARITY-THROUGH-TIME
#===========================
# 5.1 GLOBAL
#############
# Use Urtzi function
library(geiger)
source('C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/Scripts/dtt_full.R')
dttall <- dttFullCIs(tr74, svl, nsims=1000)
dttall$Pvalue
dttall$MDI


# 5.2 By radiation
###################
# We need to use the reduced dataset that does not include species from eurasia distributed in america
dt68 <- read.csv("DATA/dt_68sp.csv")
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

# EURASIA (import phylogeny of sp with svl data)
#=========
trEA <- read.tree("DATA/phyloEA_31sp.tre")
svlEA <- setNames(log(dtEA$SVL), dtEA$sp)
times <- c(50, branching.times(trEA))

dtt_EA <- dttFullCIs(trEA, svlEA, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dtt_EA$Pvalue
dtt_EA$MDI

plot(dtt_EA$dtt.data$times, dtt_EA$dtt.data$dtt, type = "l",lwd = 2, 
     xlab = "Relative time", ylab = "Disparity", axes = F, xlim=c(0,1), ylim=c(0,2))
axis(1, at= c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c("0","10", "20", "30", "40", "50"))
axis(2, at=c(0,0.5,1,1.5,2.0), labels=c("0","0.5","1.0","1.5","2.0"))
box()


#AMERICA
#========
trAM <- read.tree("DATA/phyloAM_37sp.tre")
trAM
svlAM <- setNames(log(dtAM$SVL), dtAM$sp)

dttAM <- dttFullCIs(trAM, svlAM, nsims=10000, color="#FFFFD5", linecol="gold")
dttAM$Pvalue
dttAM$MDI

plot(dttAM$times, dttAM$dtt.data$dtt, type = "l",lwd = 2, xlab = "Relative time", ylab = "Disparity")


# ALL THREE
layout(matrix(c(1,2,3), nrow = 1, ncol=3))
dttall <- dttFullCIs(tr74, svl, nsims=1000)
dtt_EA <- dttFullCIs(trEA, svlEA, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dttAM <- dttFullCIs(trAM, svlAM, nsims=10000, color="#FFFFD5", linecol="gold")

# GET TIMES
branching.times(trEA)  
dtt_EA$times
#0.6 relative time ~ 18 Mya
#0.76 relative time ~ 10 Mya

branching.times(trAM)
dttAM$times
# 0 relative time ~37 Mya
# 0.62 relative time ~ 14 Mya
# 0.9 relative time ~ 3.5 Mya

# PLOT TREES
layout(1,1)
trEA$tip.label <- gsub("Rana_", "R. ", trEA$tip.label)
plot(trEA, edge.col="darkorchid", cex=0.6)
axisPhylo(1, las = 1, lwd=2)

trAM$tip.label <- gsub("Rana_", "R. ", trAM$tip.label)
plot(trAM, edge.col="gold", cex=0.6)
axisPhylo(1, las = 1)

plot(tr, cex=0.6)
axisPhylo()


# 6. DTTs on other morphological traits
#=======================================
morpho <- read.csv("DATA/dt_morpho_55sp.csv")
morpho <- morpho[,-1]
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]
rownames(morpho) <- morpho$sp


# Separate radiations
morpho.EA <- morpho[morpho$phylo_reg=="EA",]
morpho.AM <- morpho[morpho$phylo_reg=="AM",]


# 8.1 HL
#========
library(castor)
mat.hl <- morpho[complete.cases(morpho$res_HL),]
phylo.hl <- (get_subtree_with_tips(tr74, mat.hl$sp))$subtree
HL <- setNames(mat.hl$res_HL, mat.hl$sp)

mat.hl.ea <- morpho.EA[complete.cases(morpho.EA$res_HL),]
phylo.hl.ea <- (get_subtree_with_tips(tr74, mat.hl.ea$sp))$subtree
HL.ea <- setNames(mat.hl.ea$res_HL, mat.hl.ea$sp)

mat.hl.am <- morpho.AM[complete.cases(morpho.AM$res_HL),]
phylo.hl.am <- (get_subtree_with_tips(tr74, mat.hl.am$sp))$subtree
HL.am <- setNames(mat.hl.am$res_HL, mat.hl.am$sp)


# 8.2 HW
#========
mat.hw <- morpho[complete.cases(morpho$res_HW),]
phylo.hw <- (get_subtree_with_tips(tr74, mat.hw$sp))$subtree
HW <- setNames(mat.hw$res_HW, mat.hw$sp)

mat.hw.ea <- morpho.EA[complete.cases(morpho.EA$res_HW),]
phylo.hw.ea <- (get_subtree_with_tips(tr74, mat.hw.ea$sp))$subtree
HW.ea <- setNames(mat.hw.ea$res_HW, mat.hw.ea$sp)

mat.hw.am <- morpho.AM[complete.cases(morpho.AM$res_HW),]
phylo.hw.am <- (get_subtree_with_tips(tr74, mat.hw.am$sp))$subtree
HW.am <- setNames(mat.hw.am$res_HW, mat.hw.am$sp)


# 8.3 FL
#========
mat.fl <- morpho[complete.cases(morpho$res_FL),]
phylo.fl <- (get_subtree_with_tips(tr74, mat.fl$sp))$subtree
FL <- setNames(mat.fl$res_FL, mat.fl$sp)

mat.fl.ea <- morpho.EA[complete.cases(morpho.EA$res_FL),]
phylo.fl.ea <- (get_subtree_with_tips(tr74, mat.fl.ea$sp))$subtree
FL.ea <- setNames(mat.fl.ea$res_FL, mat.fl.ea$sp)

mat.fl.am <- morpho.AM[complete.cases(morpho.AM$res_FL),]
phylo.fl.am <- (get_subtree_with_tips(tr74, mat.fl.am$sp))$subtree
FL.am <- setNames(mat.fl.am$res_FL, mat.fl.am$sp)


# 8.4 TL
#========
mat.tl <- morpho[complete.cases(morpho$res_TL),]
phylo.tl <- (get_subtree_with_tips(tr74, mat.tl$sp))$subtree
TL <- setNames(mat.tl$res_TL, mat.tl$sp)

mat.tl.ea <- morpho.EA[complete.cases(morpho.EA$res_TL),]
phylo.tl.ea <- (get_subtree_with_tips(tr74, mat.tl.ea$sp))$subtree
TL.ea <- setNames(mat.tl.ea$res_TL, mat.tl.ea$sp)

mat.tl.am <- morpho.AM[complete.cases(morpho.AM$res_TL),]
phylo.tl.am <- (get_subtree_with_tips(tr74, mat.tl.am$sp))$subtree
TL.am <- setNames(mat.tl.am$res_TL, mat.tl.am$sp)


# PLOT all together (THIS loop takes a while...)
#===================
layout(matrix(1:15, ncol=3, nrow=5))
par(mar=c(2,2,2,2))

dttall <- dttFullCIs(tr74, svl, nsims=1000)
dtt.hl <- dttFullCIs(phylo.hl, HL, nsims=10000, linecol="black")
dtt.hw <- dttFullCIs(phylo.hw, HW, nsims=10000, linecol="black")
dtt.fl <- dttFullCIs(phylo.fl, FL, nsims=10000, linecol="black")
dtt.tl <- dttFullCIs(phylo.tl, TL, nsims=10000, linecol="black")
dtt_EA <- dttFullCIs(trEA, svlEA, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dtt.hl.ea <- dttFullCIs(phylo.hl.ea, HL.ea, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dtt.hw.ea <- dttFullCIs(phylo.hw.ea, HW.ea, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dtt.fl.ea <- dttFullCIs(phylo.fl.ea, FL.ea, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dtt.tl.ea <- dttFullCIs(phylo.tl.ea, HL.ea, nsims=10000, linecol="darkorchid", color="#EBE8FC")
dttAM <- dttFullCIs(trAM, svlAM, nsims=10000, color="#FFFFD5", linecol="gold")
dtt.hl.am <- dttFullCIs(phylo.hl.am, HL.am, nsims=10000, color="#FFFFD5", linecol="gold")
dtt.hw.am <- dttFullCIs(phylo.hw.am, HW.am, nsims=10000, color="#FFFFD5", linecol="gold")
dtt.fl.am <- dttFullCIs(phylo.fl.am, FL.am, nsims=10000, color="#FFFFD5", linecol="gold")
dtt.tl.am <- dttFullCIs(phylo.tl.am, TL.am, nsims=10000, color="#FFFFD5", linecol="gold")

dtt.tl.ea$MDI
dtt.tl.ea$Pvalue


