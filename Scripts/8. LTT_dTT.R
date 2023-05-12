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
plot(ltt(phyloEA), show.tree = T, main="European radiation")
plot(ltt(phyloAM), show.tree = T, main="American radiation")

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
plot(ltt_AM$times, log(ltt_AM$ltt), type="l", col="black", lwd=2, 
     main="Lineage-through-time using all sp", xlab="time", ylab="log(nº of lineages)", xlim=c(0,50))
lines(ltt_EA$times, log(ltt_EA$ltt), type="l", col="darkgrey", lwd=2)
legend("bottomright", legend = c("American radiation", "European radiation"), 
       col = c("black", "darkgrey"), lty = 1, lwd=2, bty="n")

# No logaritm
plot(ltt_AM$times, ltt_AM$ltt, type="p", pch=16, bg="black", lwd=2, 
     main="Lineage-through-time with all sp", xlab="time", ylab="nº of lineages)", xlim=c(0,50))
lines(ltt_EA$times, ltt_EA$ltt, type="p", pch=21, bg="darkgrey", lwd=2)
legend("bottomright", legend = c("American radiation", "European radiation"), 
       pt.bg = c("black", "darkgrey"), bty="n", pch=21)


# 5. DISPARITY-THROUGH-TIME
#===========================
# 5.1 GLOBAL
#############
library(geiger)
disparity(tr74, svl)
dtt <- dtt(tr74, svl, nsim = 10000, calculateMDIp=T)
dtt$MDI
dtt$MDIpVal

# Use Urtzi function
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

dtt_EA <- dttFullCIs(trEA, svlEA, nsims=10000)
dtt_EA$Pvalue
dtt_EA$MDI

plot(dtt_EA$dtt.data$times, dtt_EA$dtt.data$dtt, type = "l",lwd = 2, 
     xlab = "Relative time", ylab = "Disparity", axes = F, xlim=c(0,1), ylim=c(0,2))
axis(1, at= c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c("0","10", "20", "30", "40", "50"))
axis(2, at=c(0,0.5,1,1.5,2.0), labels=c("0","0.5","1.0","1.5","2.0"))
box()
plot(dtt_EA$dtt.data$dtt, type="l")



# PLOT dtt togheter with LTT (but not using all species)
lttEA <- ltt(trEA)
length(lttEA$ltt)
length(dtt_EA$dtt.data$dtt)
length(dtt_EA$times)
x <- branching.times(trEA)
x <- c(x,0)
x <- sort(x)
length(x)
dtt_EA$times <- x

mat <- cbind(dtt_EA$times, dtt_EA$dtt.data$dtt)
colnames(mat) <- c("times", "dtt")

mat <- cbind(mat,x)
plot(mat[,3], mat[,2], type="l")


x <- dttFullCIs(tr, svl, nsims=10000)
plot(lttEA$times, log(lttEA$ltt), pch=21, bg="black", add=T)
lines(dtt_EA$dtt.data$dtt)

dtt_EA$times <- c(dtt_EA$times, 1)
dtt_EA$times <- lttEA$times


#AMERICA
#========
trAM <- read.tree("DATA/phyloAM_37sp.tre")
trAM
svlAM <- setNames(log(dtAM$SVL), dtAM$sp)

dttAM <- dttFullCIs(trAM, svlAM, nsims=10000)
dttAM$Pvalue
dttAM$MDI

plot(dttAM$times, dttAM$dtt.data$dtt, type = "l",lwd = 2, xlab = "Relative time", ylab = "Disparity")


# PLOT BOTH TOGHETER
plot(dtt_EA$times, dtt_EA$dtt.data$dtt, type = "l",lwd = 2, xlab = "Relative time", 
     ylim=c(0,2),ylab = "Disparity", col="darkgrey", main="Disparity-through-time in body size")
lines(dttAM$times, dttAM$dtt.data$dtt,lwd = 2)
legend("bottomleft", legend = c("American radiation", "European radiation"), 
       col = c("black", "darkgrey"), lty = 1, lwd=2, bty="n")


# 6. DTT on other morphological traits
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

# Remove sp with NAs
morpho.red <- na.omit(morpho[,c(1:4,6:11)])  # 34 sp in total

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
mat.red.EA <- mat.red[mat.red$phylo_reg=="EA",]
mat.red.AM <- mat.red[mat.red$phylo_reg=="AM",]

# Clip phylogeny
library(castor)
phylo.red.EA <- get_subtree_with_tips(trEA, mat.red.EA$sp)
phylo.red.EA <- phylo.red.EA$subtree  # 21 

plot(phylo.red.EA)
axisPhylo()  # 50 MYa

phylo.red.AM <- get_subtree_with_tips(trAM, mat.red.AM$sp)
phylo.red.AM <- phylo.red.AM$subtree  # 13 species

plot(phylo.red.AM)
axisPhylo() # 30 Mya


# 8.1 HL
#========
HL.ea <- setNames(mat.red.EA$res_HL, mat.red.EA$sp)
dtt.hl.ea <- dttFullCIs(phylo.red.EA, HL.ea, nsims=10000, linecol="darkorchid")
dtt.hl.ea$Pvalue
dtt.hl.ea$MDI

HL.am <- setNames(mat.red.AM$res_HL, mat.red.AM$sp)
dtt.hl.am <- dttFullCIs(phylo.red.AM, HL.am, nsims=10000, linecol="gold")
dtt.hl.am$Pvalue
dtt.hl.am$MDI

# 8.2 HW
#========
HW.ea <- setNames(mat.red.EA$res_HW, mat.red.EA$sp)
dtt.hw.ea <- dttFullCIs(phylo.red.EA, HW.ea, nsims=10000, linecol="darkorchid")
dtt.hw.ea$Pvalue
dtt.hw.ea$MDI

HW.am <- setNames(mat.red.AM$res_HW, mat.red.AM$sp)
dtt.hw.am <- dttFullCIs(phylo.red.AM, HW.am, nsims=10000, linecol = "gold")
dtt.hw.am$Pvalue
dtt.hw.am$MDI

# 8.3 FL
#========
FL.ea <- setNames(mat.red.EA$res_FL, mat.red.EA$sp)
dtt.fl.ea <- dttFullCIs(phylo.red.EA, FL.ea, nsims=10000, linecol="darkorchid")
dtt.fl.ea$Pvalue
dtt.fl.ea$MDI

FL.am <- setNames(mat.red.AM$res_FL, mat.red.AM$sp)
dtt.fl.am <- dttFullCIs(phylo.red.AM, FL.am, nsims=10000, linecol = "gold")
dtt.fl.am$Pvalue
dtt.fl.am$MDI

# 8.4 TL
#========
TL.ea <- setNames(mat.red.EA$res_TL, mat.red.EA$sp)
dtt.tl.ea <- dttFullCIs(phylo.red.EA, TL.ea, nsims=10000, linecol="darkorchid")
dtt.tl.ea$Pvalue
dtt.tl.ea$MDI

TL.am <- setNames(mat.red.AM$res_TL, mat.red.AM$sp)
dtt.tl.am <- dttFullCIs(phylo.red.AM, TL.am, nsims=10000, linecol = "gold")
dtt.tl.am$Pvalue
dtt.tl.am$MDI

plot(dtt.tl.am$times, dtt.tl.am$dtt.data$dtt, type="l")


# 9. All morphological variables togheter
#=========================================
mat.red.AM <- mat.red.AM[phylo.red.AM$tip.label,]
mat.red.AM$sp == phylo.red.AM$tip.label

mat.red.EA <- mat.red.EA[phylo.red.EA$tip.label,]
mat.red.EA$sp == phylo.red.EA$tip.label

morpho.AM <- mat.red.AM[,11:14]
morpho.EA <- mat.red.EA[,11:14]

# EURASIA
dtt.morpho.ea <- dttFullCIs(phylo.red.EA, morpho.EA, nsims=10000, linecol = "darkorchid")
dtt.morpho.ea$MDI
dtt.morpho.ea$Pvalue

# AMERICA
dtt.morpho.am <- dttFullCIs(phylo.red.AM, morpho.AM, nsims=10000, linecol = "gold")
dtt.morpho.am$MDI
dtt.morpho.am$Pvalue



####################################################################
# 10.ML methods for analyzing time shifts in diversification rates  
####################################################################
library(laser)
?fitdAICrc  # Test for Rate Variation Using delta-AICrc Test Statistic
?DDX        # Fit Density Dependent Speciation Model to Branching Times
?yule2rate  # Fits multi-rate variants of the pure birth (Yule) model 
# in the case of yule2rate, allows for one shift into a new rate. 

# GLOBAL 
#=========
alltimes <- branching.times(tr)
alltimes
div.models <- fitdAICrc(alltimes, modelset = c("pureBirth", "bd","DDX", "DDL", "yule2rate", "yule3rate"), ints = 100) 

# EURASIA
#==========
EAtimes <- branching.times(phyloEA)
EAtimes
div.models.EA <- fitdAICrc(EAtimes, modelset = c("pureBirth", "bd","DDX", "DDL", "yule2rate", "yule3rate"), ints = 100) 
# Best model with the lowest AIC = 70 --> yule2rate model


# AMERICA
#=========
AMtimes <- branching.times(phyloAM)
AMtimes
div.models.AM <- fitdAICrc(AMtimes, modelset = c("pureBirth", "bd","DDX", "DDL", "yule2rate", "yule3rate"), ints = 100) 
# Best model with the lowest AIC = 70 --> yule2rate model


# SAVE
# write.csv(div.models, "Results/div.models.csv")
# write.csv(div.models.AM, "Results/div.models.AM.csv")
# write.csv(div.models.EA, "Results/div.models.EA.csv")
