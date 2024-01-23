########################################################################
# STRAPP: Structured rate permutations on phylogenies 
#  Assessing correlation between a continuous trait and speciation rate
#   using permutations
########################################################################
# For this analysis we need the output of BAMM

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(phytools)
library(ape)
library(BAMMtools)


#################
# GLOBAL (68 sp)
#################

# 1. Import data
#================
tr74 <- read.tree('DATA/phyloRana_74sp.tre')
tr68 <- read.tree("DATA/phyloRana_68sp.tre")

# BAMM output diversification
ediv <- getEventData(tr68, eventdata = "BAMM/bamm-2.5.0-Windows/event_data_div68.txt", 
                     burnin=0.1, type = "diversification")

# BAMM output phenotypic evolution
etrait <- getEventData(tr68, eventdata = "BAMM/bamm-2.5.0-Windows/event_data_trait_68.txt", 
                               burnin=0.5, type = "trait")

# 2. Assign rates to each tip of the phylogeny
#=============================================
# rates of phenotypic evolution, speciation and net diversification 
traits <- getTipRates(etrait) 
div <- getTipRates(ediv, returnNetDiv = T)
spec <- getTipRates(ediv, returnNetDiv = F)

traits <- traits$beta.avg  # Average tip phenotypic rates
div <- div$netdiv.avg      # Average tip diversification rates
spec <- spec$lambda.avg    # Average tip speciation rates 

#3. Permutation analysis
#=======================
# Diversification ~ phenotypic rate
rana.permu.sp<-traitDependentBAMM(ediv, traits, reps=10000, return.full = T,  
                               two.tailed = T, method='s')
rana.permu.sp$estimate # Very low correlation between speciation rates and svl evolution
rana.permu.sp$p.value  # NOT significant

par(mar=c(6,6,6,6))
plot(y=div, x=traits, pch=21, bg="black", ylab="Net diversification rate", xlab="phenotypic rate (SVL)")
abline(lm(div~traits), col="red", lwd=2)
text(x=0.03, y=0.05, "Global (68sp)", font=4)

cor(traits, div)


#################
# EURASIA (31 sp)
#################
# 1. Import data
#================
trEA <- read.tree('DATA/phyloEA_31sp.tre')

# BAMM output diversification
edivEA <- getEventData(trEA, eventdata = "BAMM/bamm-2.5.0-Windows/OLD/event_data_div_EA.txt", burnin=0.1)

# BAMM output phenotypic evolution
etraitEA <- getEventData(trEA, eventdata = "BAMM/bamm-2.5.0-Windows/OLD/event_data_trait_EA_200.txt", 
                       burnin=0.5, type = "trait")

# 2. Assign rates to each tip of the phylogeny
#=============================================
# rates of sp, extinction and net diversification 
traitsEA <- getTipRates(etraitEA)$beta.avg  # Average tip phenotypic rates 
divEA <- getTipRates(edivEA, returnNetDiv = T)$netdiv.avg     # Average tip diversification rates
specEA <- getTipRates(edivEA)$lambda.avg    # Average tip speciation rates

# Diversification ~ phenotypic rate
ea.permu<-traitDependentBAMM(edivEA, traitsEA, 10000, return.full = T, logrates = T, 
                               two.tailed = T, method='s')
ea.permu$estimate # Negative correlation between diversification rates and svl evolution
ea.permu$p.value  # NOT significant

par(mar=c(6,6,6,6))
plot(x=traitsEA, y=divEA, pch=21, bg="black", ylab="Diversification rate", xlab="Phenotypic rate (SVL)")
abline(lm(divEA~traitsEA), col="red", lwd=2)
#text(x=0.0062, y=0, "Eurasia (31sp)", font=4, cex=1.3)


#################
# AMERICA (37 sp)
#################
# 1. Import data
#================
trAM <- read.tree('DATA/phyloAM_37sp.tre')

# BAMM output diversification
edivAM <- getEventData(trAM, eventdata = "BAMM/bamm-2.5.0-Windows/OLD/event_data_div_AM.txt", burnin=0.1)

# BAMM output phenotypic evolution
etraitAM <- getEventData(trAM, eventdata = "BAMM/bamm-2.5.0-Windows/OLD/event_data_trait_AM_200.txt", 
                         burnin=0.5, type = "trait")

# 2. Assign rates to each tip of the phylogeny
#=============================================
# rates of sp, extinction and net diversification 
traitsAM <- getTipRates(etraitAM)$beta.avg  # Average tip phenotypic rates 
divAM <- getTipRates(edivAM, returnNetDiv = T)$netdiv.avg     # Average tip diversification rates
specAM <- getTipRates(edivAM)$lambda.avg    # Average tip speciation rates

# Diversification ~ phenotypic rate
am.permu<-traitDependentBAMM(edivAM, traitsAM, 10000, return.full = T, logrates = T, 
                             two.tailed = T, method='s')
am.permu$estimate # Negative correlation between diversification rates and svl evolution
am.permu$p.value  # NOT significant

par(mar=c(6,6,6,6))
plot(x=traitsAM, y=divAM, pch=21, bg="black", ylab="Diversification rate", xlab="Phenotypic rate (SVL)")
abline(lm(divAM~traitsAM), col="red", lwd=2)
# text(x=0.04, y=0, "America (37sp)", font=4, cex=1.3)




#PLOT 3 (GGPLOT)
library(ggplot2)
mat.all <- as.data.frame(cbind(traits, div))
mat.ea <- as.data.frame(cbind(traitsEA, divEA))
mat.am <- as.data.frame(cbind(traitsAM, divAM))

pa <- ggplot(mat.all, aes(x=traits, y=div)) + 
  geom_point() +
  theme_classic() +
  ylab("Diversification Rate") +
  xlab("Rate of SVL evolution") +
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Global")

pb <- ggplot(mat.ea, aes(x=traitsEA, y=divEA)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of SVL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Eurasia")

pc <- ggplot(mat.am, aes(x=traitsAM, y=divAM)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of SVL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") +
  ggtitle("America")

library(patchwork)
pa + pb + pc


