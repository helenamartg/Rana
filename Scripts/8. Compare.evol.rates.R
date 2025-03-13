###############################
#     COMPARE EVOL RATES
# of size-corrected variables
###############################

rm(list=ls())

library(phytools)
library(RRPP)
library(geiger)
library(pls)
library(geomorph)
library(castor)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")


# 1. Read phylogeny, morphometrics and environmental variables, and separate genera
#==================================================================================
phylorana_45 <- read.tree("Data/phyloRana_45sp.tre")

traits <- read.csv("Data/morpho_phylores_45sp.csv")
rownames(traits) <- traits$Species
morpho <- traits[,c(2,4,12:16)]

morpho <- na.omit(morpho)
phylo28 <- (get_subtree_with_tips(phylorana_45, morpho$Species))$subtree



# 2. Check that phenotypic data and tree tips are in the same order
#==================================================================
phylo28$tip.label ==  morpho$Species
morpho <- morpho[phylo28$tip.label,]


# 3. Compare evol.rates
#=======================
rad <- setNames(morpho$phylo_reg, morpho$Species)
morph2 <- morpho[,4:7]

evolrate.morpho <- compare.evol.rates(morph2, phylo28, rad, method = "simulation")
evolrate.morpho

evolrate.morpho$sigma.d.gp


plot(evolrate.morpho)
plot(evolrate.morpho$sigma.d.gp)
evo <- cbind(evolrate.morpho$groups, round(evolrate.morpho$sigma.d.gp, 5))
colnames(evo) <- c("groups", "sigma")
evo <- as.data.frame(evo)

# plot
library(ggplot2)
ggplot(evo, aes(x=groups, y=sigma)) +
  theme_classic() +
  ylab("sigma") +
  geom_point(size=5, colour=c("gold", "darkorchid"))



# Compare evol rates of body size
#==================================
dt59 <- read.csv("DATA/dt_59sp.csv")
phylo59 <- read.tree("DATA/phyloRana_59sp.tre")

rownames(dt59) <- dt59$Species
dt59 <- dt59[phylo59$tip.label,]
svl <- setNames(dt59$log_svl, dt59$Species)
evolrate.svl <- compare.evol.rates(svl, phylo59, rad, method="simulation")

evolrate.svl
round(evolrate.svl$sigma.d.gp, 3)

plot(evolrate.svl)
plot(evolrate.svl$sigma.d.gp)
evo.svl <- cbind(evolrate.svl$groups, round(evolrate.svl$sigma.d.gp, 5))
colnames(evo.svl) <- c("groups", "sigma")
evo.svl <- as.data.frame(evo.svl)

ggplot(evo.svl, aes(x=groups, y=sigma)) +
  theme_classic() +
  ylab("sigma") +
  geom_point(size=5, colour=c("gold", "darkorchid"))


