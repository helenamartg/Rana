########################################################################
# STRAPP: Structured rate permutations on phylogenies 
#  Assessing correlation between a continuous trait and speciation rate
#   using permutations
########################################################################
# For this analysis we need the output of BAMM

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

library(phytools)
library(ape)
library(BAMMtools)


#################
# GLOBAL (65 sp)
#################

# 1. Import data
#================
tr59 <- read.tree("Data/phyloRana_59sp.tre")
tr.svl<- read.tree('BAMM/tree_svl.tre')

# BAMM output diversification
ediv <- getEventData(tr59, eventdata = "BAMM/event_data_div59.txt", 
                     burnin=0.1, type = "diversification")

# BAMM output phenotypic evolution
etrait <- getEventData(tr.svl, eventdata = "BAMM/event_data_trait_svl.txt", 
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
                               two.tailed = T, method='s', logrates = T)
rana.permu.sp$estimate # Very low correlation between speciation rates and svl evolution
rana.permu.sp$p.value  # NOT significant

par(mar=c(6,6,6,6))
plot(y=div, x=traits, pch=21, bg="black", ylab="Net diversification rate", xlab="phenotypic rate (SVL)")
abline(lm(div~traits), col="red", lwd=2)
text(x=0.03, y=0.05, "Global (68sp)", font=4)

cor(traits, div)


### SEPARATE RADIATIONS ###
# Extract from global
names(traits) == names(div)
global <- cbind(traits, div)
global <- as.data.frame(global)
global$sp <- rownames(global)

sp59 <- read.csv("Data/dt_59sp.csv")
rownames(sp59) <- sp59$Species

global$sp == sp59$Species
global <- global[sp59$Species,]
global$radiation <- sp59$phylo_subtree



#################
# EURASIA (26 sp)
#################
trEA <- read.tree('Data/phyloEA_26sp.tre')

rates.EA <- global[global$radiation=="EA",]

# Assign rates to each tip of the phylogeny
divEA <- setNames(rates.EA$div, rates.EA$sp)
traitsEA <- setNames(rates.EA$traits, rates.EA$sp)

edivEA <- getEventData(trEA, eventdata = "BAMM/event_data_div_EA.txt", 
                     burnin=0.1, type = "diversification")


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
# AMERICA (33 sp)
#################
trAM <- read.tree('Data/phyloAM_33sp.tre')

rates.AM <- global[global$radiation=="AM",]
rates.AM <- as.data.frame(rates.AM)

# Assign rates to each tip of the phylogeny
div.AM <- setNames(rates.AM$div, rates.AM$sp)
traitsAM <- setNames(rates.AM$traits, rates.AM$sp)

edivAM <- getEventData(trAM, eventdata = "BAMM/event_data_div_AM.txt", 
                       burnin=0.1, type = "diversification")


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


