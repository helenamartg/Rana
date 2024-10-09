#######################################
# Getting ready for BAMM analyses
######################################

rm(list=ls())

library(phytools)
library(ape)
library(castor)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/BAMM")

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
tr74 <- read.tree("bamm-2.5.0-Windows/phyloRana_74sp.tre")
is.ultrametric(tr74)   # Time-calibrated
is.binary(tr74)        # No polytomies
min(tr74$edge.length)  # Branch lengths > 0 

tr68 <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/phyloRana_68sp.tre")
is.ultrametric(tr68)



# 2. Set priors of evol. rate parameters (modeltype = speciationextinction)
#========================================
library(BAMMtools)
setBAMMpriors(tr83, outfile = "mupriors_div_83sp.txt")
setBAMMpriors(tr68, outfile = "myPriors_div_68sp.txt")


# 3. Set priors of evol. rate parameters (modeltype = trait)
#=======================================
tr68 <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/phyloRana_68sp.tre")
traits <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/dt_68sp.csv")
svl <- cbind(traits$sp, log(traits$SVL))
rownames(svl) <- svl[,1]
svl[,1]==tr68$tip.label
svl <- svl[tr68$tip.label,]
svl <- svl[,-1]
svl <- as.data.frame(svl)
# write.table(svl, "svl.txt")

setBAMMpriors(phy=tr68, traits = "svl.txt", outfile = "myPriors_traits.txt")


# 4. Set priors for the rest of morphological variables
#=======================================================
morpho <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/dt_phylores_55sp.csv")
tr55 <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/phyloRana_55sp.tre")

# Remove 6 species + R. tlaloci
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]


tr55 <- get_subtree_with_tips(tr55, morpho$sp)
tr51 <- tr55$subtree
plot(tr51)


############
# HW (51sp)
############
hw <- morpho[,c(2,15)]
rownames(hw) <- hw[,1]
hw[,1] == tr51$tip.label
hw <- hw[tr51$tip.label,]
hw <- setNames(hw$phylores_HW, hw$sp)
hw <- na.omit(hw)
hw <- as.data.frame(hw)
# write.table(hw, "hw.txt", sep=" ")
# write.tree(tr51, "tree_hw.tre")

setBAMMpriors(phy = tr51, traits = "hw.txt", outfile = "myPriors_hw.txt")


############
# HL (45sp)
############
hl <- morpho[,c(2,14)]
hl[,1] == tr51$tip.label
rownames(hl) <- hl$sp
hl <- as.data.frame(hl)
hl <- setNames(hl$phylores_HL, hl$sp)
hl <- na.omit(hl)
# write.table(hl, "hl.txt")

tr51
tr45 <- get_subtree_with_tips(tr51, names(hl))
tr45 <- tr45$subtree
# write.tree(tr45, "tree_hl.tre")
names(hl) == tr45$tip.label  # CHECK

# Get priors
setBAMMpriors(phy = tr45, traits = "hl.txt", outfile = "myPriors_hl.txt")



############
# FL (43sp)
############
fl <- morpho[,c(2,16)]
rownames(fl) <- fl[,1]
fl[,1] == tr51$tip.label
fl <- fl[tr51$tip.label,]
fl <- setNames(fl$phylores_FL, fl$sp)
fl <- na.omit(fl)
# write.table(fl, "fl.txt")

tr51
tr43 <- get_subtree_with_tips(tr51, names(fl))
tr43 <- tr43$subtree
# write.tree(tr43, "tree_fl.tre")
names(fl) == tr43$tip.label  # CHECK

# Get priors
setBAMMpriors(phy = tr43, traits = "fl.txt", outfile = "myPriors_fl.txt")



############
# TL (46sp)
############
tl <- morpho[,c(2,17)]
rownames(tl) <- tl[,1]
tl[,1] == tr51$tip.label
tl <- tl[tr51$tip.label,]
tl <- setNames(tl$phylores_TL, tl$sp)
tl <- na.omit(tl)
# write.table(tl, "tl.txt")

tr51
tr46 <- get_subtree_with_tips(tr51, names(tl))
tr46 <- tr46$subtree
# write.tree(tr46, "tree_tl.tre")
names(tl) == tr46$tip.label  # CHECK

# Get priors
setBAMMpriors(phy = tr46, traits = "tl.txt", outfile = "myPriors_tl.txt")

