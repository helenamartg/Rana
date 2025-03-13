#######################################
# Getting ready for BAMM analyses
######################################

rm(list=ls())

library(phytools)
library(ape)
library(castor)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4publication/2nd_Revision/Reanalisis")

# 1. Import the tree and check BAMM assumptions
#===============================================
# YUAN phylogeny
tr <- read.nexus("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/raw_data/Yuan_2016_phylo_v2.tree")
is.ultrametric(tr)
par(mar=c(1,1,1,1))
plot(tr, cex=0.6)
nodelabels()

tr <- extract.clade(tr, 97) # Remove outgroups
plot(tr, cex=0.6) #New tree has 90 tips

tr$tip.label # Remove R_kunyuensis (= R. coreana)
tr <- drop.tip(tr, 2) #New tree with 89 tips
plot(tr, cex=.6)
is.ultrametric(tr)
tr83 <- drop.clade(tr, c(141:144))
tr83 <- drop.tip(tr83, "NA")
plot(tr83, cex=.6)
# write.tree(tr83, "phyloRana_83sp.tre")

# Import trees without undescribed sp
tr59 <- read.tree("Data/phyloRana_59sp.tre")
is.ultrametric(tr59)   # Time-calibrated
is.binary(tr59)        # No polytomies
min(tr59$edge.length)  # Branch lengths > 0 

trEA <- read.tree("Data/phyloEA_26sp.tre")
trAM <- read.tree("Data/phyloAM_33sp.tre")

is.ultrametric(trEA)
is.ultrametric(trAM)


# 2. Set priors of evol. rate parameters (modeltype = speciationextinction)
#========================================
library(BAMMtools)
# setBAMMpriors(tr83, outfile = "BAMM&mypriors_div_83sp.txt")
setBAMMpriors(tr59, outfile = "BAMM/myPriors_div_59sp.txt")
setBAMMpriors(trEA, outfile="BAMM/mypriors_div_EA.txt")
setBAMMpriors(trAM, outfile="BAMM/mypriors_div_AM.txt")



# 3. Set priors of evol. rate parameters (modeltype = trait)
#=======================================
traits <- read.csv("Data/dt_59sp.csv")
rownames(traits) <- traits$Species
svl <- traits[,c(2,9)]

# check order
svl$Species==tr59$tip.label
svl <- svl[tr59$tip.label,]
svl <- svl[,-1]
svl <- as.data.frame(svl)
rownames(svl) <- traits$Species

# write.table(svl, "BAMM/svl.txt")
# setBAMMpriors(phy=tr59, traits = "BAMM/svl.txt", outfile = "BAMM/myPriors_svl.txt")


# 4. Set priors for the rest of morphological variables
#=======================================================
morpho <- read.csv("Data/morpho_phylores_45sp.csv")
tr45 <- read.tree("Data/phyloRana_45sp.tre")

# Remove 6 species + R. tlaloci
morpho <- morpho[morpho$Species!="Rana_tlaloci",]
morpho <- morpho[morpho$Species!="Rana_luteiventris",]
morpho <- morpho[morpho$Species!="Rana_boylii",]
morpho <- morpho[morpho$Species!="Rana_sierrae",]
morpho <- morpho[morpho$Species!="Rana_muscosa",]
morpho <- morpho[morpho$Species!="Rana_aurora",]
morpho <- morpho[morpho$Species!="Rana_cascadae",]


tr42 <- get_subtree_with_tips(tr45, morpho$Species)$subtree
plot(tr42)
is.ultrametric(tr42)


############
# HW (42sp)
############
hw <- morpho[,c(2,14)]
rownames(hw) <- hw[,1]
hw[,1] == tr42$tip.label
hw <- hw[tr42$tip.label,]
hw <- setNames(hw$phylores_HW, hw$Species)
hw <- na.omit(hw)
hw <- as.data.frame(hw)
write.table(hw, "BAMM/hw.txt", sep="\t")
write.tree(tr42, "BAMM/tree_hw.tre")

setBAMMpriors(phy = tr42, traits = "BAMM/hw.txt", outfile = "BAMM/myPriors_hw.txt")


############
# HL (36sp)
############
hl <- morpho[,c(2,13)]
rownames(hl) <- hl[,1]
hl[,1] == tr42$tip.label
hl <- hl[tr42$tip.label,]
hl <- setNames(hl$phylores_HL, hl$Species)
hl <- na.omit(hl)
hl <- as.data.frame(hl)
# write.table(hl, "BAMM/hl.txt", sep="\t")

tr42
tr36 <- get_subtree_with_tips(tr42, rownames(hl))$subtree
# write.tree(tr36, "BAMM/tree_hl.tre")

rownames(hl) == tr36$tip.label  # CHECK


# Get priors
# setBAMMpriors(phy = tr36, traits = "BAMM/hl.txt", outfile = "BAMM/myPriors_hl.txt")



############
# FL (35sp)
############
fl <- morpho[,c(2,15)]
rownames(fl) <- fl[,1]
fl[,1] == tr42$tip.label
fl <- fl[tr42$tip.label,]
fl <- setNames(fl$phylores_FL, fl$Species)
fl <- na.omit(fl)
fl <- as.data.frame(fl)
# write.table(fl, "BAMM/fl.txt", sep="\t")

tr35 <- get_subtree_with_tips(tr42, rownames(fl))$subtree

# write.tree(tr35, "BAMM/tree_fl.tre")
rownames(fl) == tr35$tip.label  # CHECK

# Get priors
# setBAMMpriors(phy = tr35, traits = "BAMM/fl.txt", outfile = "BAMM/myPriors_fl.txt")



############
# TL (46sp)
############
tl <- morpho[,c(2,16)]
rownames(tl) <- tl[,1]
tl[,1] == tr42$tip.label
tl <- tl[tr42$tip.label,]
tl <- setNames(tl$phylores_TL, tl$Species)
tl <- na.omit(tl)
tl <- as.data.frame(tl)
# write.table(tl, "BAMM/tl.txt", sep="\t")

tr38 <- get_subtree_with_tips(tr42, rownames(tl))$subtree
# write.tree(tr38, "BAMM/tree_tl.tre")
rownames(tl) == tr38$tip.label  # CHECK

# Get priors
# setBAMMpriors(phy = tr38, traits = "BAMM/tl.txt", outfile = "BAMM/myPriors_tl.txt")

