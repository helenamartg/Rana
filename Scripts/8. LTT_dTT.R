###############################
#   LTT and dTT
###############################

rm(list=ls())

library(phytools)
library(ape)
library(geomorph)

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

dtt_data <- data.frame(time = dtt$time, disparity = dtt$disparity,
                       lower_ci = dtt$lower, upper_ci = dtt$upper)

# Use Urtzi function
source('C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/dtt_full.R')
dttall <- dttFullCIs(tr74, svl, nsims=1000)
dttall$Pvalue


# 5.2 By radiation
###################
# We need to use the reduced dataset that does not include species from eurasia distributed in america
dt68 <- read.csv("DATA/dt_68sp.csv")
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

# EURASIA (import phylogeny of sp with svl data)
trEA <- read.tree("DATA/phyloEA_31sp.tre")
svlEA <- setNames(dtEA$SVL, dtEA$sp)
times <- c(50, branching.times(trEA))

dtt_EA <- dttFullCIs(trEA, svlEA, nsims=10000)
dtt_EA$Pvalue
dtt_EA$MDI

plot(dtt_EA$times, dtt_EA$dtt.data$dtt, type = "l",lwd = 2, 
     xlab = "Relative time", ylab = "Disparity", axes = F, xlim=c(0,1), ylim=c(0,2))
axis(1, at= c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c("0","10", "20", "30", "40", "50"))
axis(2, at=c(0,0.5,1,1.5,2.0), labels=c("0","0.5","1.0","1.5","2.0"))
box()
plot(dtt_EA$dtt.data$dtt, type="l")

#AMERICA
trAM <- read.tree("DATA/phyloAM_37sp.tre")
trAM
svlAM <- setNames(dtAM$SVL, dtAM$sp)

dttAM <- dttFullCIs(trAM, svlAM, nsims=10000)
dttAM$Pvalue
dttAM$MDI

plot(dttAM$times, dttAM$dtt.data$dtt, type = "l",lwd = 2, xlab = "Relative time", ylab = "Disparity")


# PLOT BOTH TOGHETER
plot(dttEA$times, dttEA$dtt.data$dtt, type = "l",lwd = 2, xlab = "Relative time", 
     ylim=c(0,2),ylab = "Disparity", col="darkgrey", main="Disparity-through-time in body size")
lines(dttAM$times, dttAM$dtt.data$dtt,lwd = 2)
legend("bottomleft", legend = c("American radiation", "European radiation"), 
       col = c("black", "darkgrey"), lty = 1, lwd=2, bty="n")




# MEDUSA
medusa.rana <- medusa(phyloRana)
plot(medusa.rana)

medusa.EA <- medusa(phyloEA)
plot(medusa.EA)

medusa.AM <- medusa(phyloAM)
plot(medusa.AM)




# 6.ML methods for analyzing time shifts in diversification rates  
#################################################################
library(laser)
?fitdAICrc  # Test for Rate Variation Using delta-AICrc Test Statistic
?DDX        # Fit Density Dependent Speciation Model to Branching Times
?yule2rate  # Fits multi-rate variants of the pure birth (Yule) model 
# in the case of yule2rate, allows for one shift into a new rate. 

EAtimes <- branching.times(phyloEA)
EAtimes
div.models.EA <- fitdAICrc(EAtimes, modelset = c("pureBirth", "bd","DDX", "DDL", "yule2rate", "yule3rate"), ints = 100) 

# Best model with the lowest AIC = 64.65019 --> DDL (logistic variant of the density-dependent speciation rate model)

