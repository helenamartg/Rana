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

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")


# 1. Read phylogeny, morphometrics and environmental variables, and separate genera
#==================================================================================
phylorana_55 <- read.tree("DATA/phyloRana_55sp.tre")

traits <- read.csv("DATA/dt_phylores_55sp.csv")
rownames(traits) <- traits$sp
morpho <- traits[,c(2,4,13:17)]

morpho <- na.omit(morpho)
phylo37 <- get_subtree_with_tips(phylorana_55, morpho$sp)
phylo37 <- phylo37$subtree


# 2. Check that phenotypic data and tree tips are in the same order
#==================================================================
phylo37$tip.label ==  morpho$sp
morpho <- morpho[phylo37$tip.label,]


# 3. Compare evol.rates
#=======================
rad <- setNames(morpho$phylo_reg, morpho$sp)
morph2 <- morpho[,4:7]

evolrate.morpho <- compare.evol.rates(morph2, phylo37, rad, method = "simulation")
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
dt68 <- read.csv("DATA/dt_68sp.csv")
phylo68 <- read.tree("DATA/phyloRana_68sp.tre")

rownames(dt68) <- dt68$sp
dt68 <- dt68[phylo68$tip.label,]
svl <- setNames(dt68$log_svl, dt68$sp)
evolrate.svl <- compare.evol.rates(svl, phylo68, rad, method="simulation")

evolrate.svl
evolrate.svl $sigma.d.gp

plot(evolrate.svl)
plot(evolrate.svl$sigma.d.gp)
evo.svl <- cbind(evolrate.svl$groups, round(evolrate.svl$sigma.d.gp, 5))
colnames(evo.svl) <- c("groups", "sigma")
evo.svl <- as.data.frame(evo.svl)

ggplot(evo.svl, aes(x=groups, y=sigma)) +
  theme_classic() +
  ylab("sigma") +
  geom_point(size=5, colour=c("gold", "darkorchid"))


