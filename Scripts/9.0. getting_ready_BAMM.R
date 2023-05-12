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
tr74 <- read.tree("bamm-2.5.0-Windows/phyloRana_74sp.tre")

is.ultrametric(tr74)   # Time-calibrated
is.binary(tr74)        # No polytomies
min(tr74$edge.length)  # Branch lengths > 0 

# Separate radiations
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

phyloAM <- extract.clade(tr, 92)
plot(phyloAM, cex=0.7)  # 47 sp (8 undescribed species)

tr$tip.label
phyloEA <- drop.clade(tr, c(93:137))
plot(phyloEA, cex=0.7)
phyloEA$tip.label
phyloEA <- drop.tip(phyloEA, "NA")
plot(phyloEA)    # 42 sp
nodelabels()

# Remove 6 species from phyloEA distributed in AM (node 46)
phyloEA <- drop.clade(phyloEA, c(47:50))
phyloEA <- drop.tip(phyloEA, "NA")
plot(phyloEA)  # 36 sp (42-6)

is.ultrametric(phyloAM)
is.ultrametric(phyloEA)

# SAVE
# write.tree(tr, "Yuan_89sp.tre")
# write.tree(phyloAM, "phyloAM_47sp.tre")
# write.tree(phyloEA, "phyloEA_36sp.tre")


# Import trees without undescribed sp
phyloEA_red <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/phyloEA_31sp.tre")
phyloAM_red <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/phyloAM_37sp.tre")

is.ultrametric(phyloAM_red)
is.ultrametric(phyloEA_red)


# 2. Set priors of evol. rate parameters (modeltype = speciationextinction)
#========================================
library(BAMMtools)
# setBAMMpriors(tr, outfile = "myPriors_diversification_ALL89sp.txt")
setBAMMpriors(phyloEA_red, outfile = "myPriors_diversification_EAred.txt")
setBAMMpriors(phyloAM_red, outfile = "myPriors_diversification_AMred.txt")
# setBAMMpriors(phyloEA, outfile = "myPriors_diversification_EA.txt")
# setBAMMpriors(phyloAM, outfile = "myPriors_diversification_AM.txt")


# 3. Import traits
#===================
traits <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/dt_74sp.csv")
svl <- cbind(traits$sp, log(traits$SVL))
rownames(svl) <- svl[,1]
svl[,1]==tr$tip.label
svl <- svl[tr$tip.label,]
svl <- svl[,-1]
svl <- as.data.frame(svl)
# write.table(svl, "svl.txt")


# 4. Set priors of evol. rate parameters (modeltype = trait)
#=======================================
setBAMMpriors(phy=tr, traits = "svl.txt", outfile = "myPriors_traits.txt")


# 5. By radiation
#=================
dt68 <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/dt_68sp.csv")
trEA <- read.tree("bamm-2.5.0-Windows/phyloEA_31sp.tre")
trAM <- read.tree("bamm-2.5.0-Windows/phyloAM_37sp.tre")
  
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

svlEA <- cbind(dtEA$sp, log(dtEA$SVL))
rownames(svlEA) <- svlEA[,1]
svlEA[,1]==trEA$tip.label
svlEA<- svlEA[trEA$tip.label,]
svlEA <- svlEA[,-1]
svlEA <- as.data.frame(svlEA)
# write.table(svlEA, "svlEA.txt")

svlAM <- cbind(dtAM$sp, log(dtAM$SVL))
rownames(svlAM) <- svlAM[,1]
svlAM[,1]==trAM$tip.label
svlAM<- svlAM[trAM$tip.label,]
svlAM <- svlAM[,-1]
svlAM <- as.data.frame(svlAM)
# write.table(svlAM, "svlAM.txt")

# Set parameters
# setBAMMpriors(phy = trEA, traits = "svlEA.txt", outfile = "myPriors_traits_EA.txt")
# setBAMMpriors(phy = trAM, traits = "svlAM.txt", outfile = "myPriors_traits_AM.txt")


# 6. Set priors for the rest of morphological variables
#=======================================================
morpho <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/dt_morpho_55sp.csv")
tr55 <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/phyloRana_55sp.tre")

#############
# HW 
#############
hw <- cbind(morpho$sp, log(morpho$HW))
rownames(hw) <- hw[,1]
hw[,1] == tr55$tip.label
hw <- hw[tr55$tip.label,]
hw <- hw[,-1]
hw <- as.data.frame(hw)
# write.table(hw, "hw.txt")

# EURASIA
plot(tr55)
nodelabels()

morphoEA <- morpho[morpho$phylo_reg=="EA",]
# Extract species from EA radiation that are distributed in aM
morphoEA <- morphoEA[morphoEA$sp!="Rana_boylii",]
morphoEA <- morphoEA[morphoEA$sp!="Rana_sierrae",]
morphoEA <- morphoEA[morphoEA$sp!="Rana_muscosa",]

phyEA <- get_subtree_with_tips(tr55, morphoEA$sp)
phyEA <- phyEA$subtree # 24 sp

hwEA <- cbind(morphoEA$sp, log(morphoEA$HW))
rownames(hwEA) <- hwEA[,1]
hwEA[,1] == phyEA$tip.label
hwEA <- hwEA[phyEA$tip.label,]
hwEA <- hwEA[,-1]
hwEA <- as.data.frame(hwEA)
# write.table(hwEA,"hwEA.txt")

# AMERICA
morphoAM <- morpho[morpho$phylo_reg=="AM",]
phyAM <- get_subtree_with_tips(tr55, morphoAM$sp)
phyAM <- phyAM$subtree # 27 sp

hwAM <- cbind(morphoAM$sp, log(morphoAM$HW))
rownames(hwAM) <- hwAM[,1]
hwAM[,1] == phyAM$tip.label
hwAM <- hwAM[phyAM$tip.label,]
hwAM <- hwAM[,-1]
hwAM <- as.data.frame(hwAM)
# write.table(hwAM,"hwAM.txt")

# SET PRIORS
setBAMMpriors(phy = tr55, traits = "hw.txt", outfile = "myPriors_hw.txt")
setBAMMpriors(phy = phyEA, traits = "hwEA.txt", outfile = "myPriors_hw_EA.txt")
setBAMMpriors(phy = phyAM, traits = "hwAM.txt", outfile = "myPriors_hw_AM.txt")

# SAVE TREES
# write.tree(phyAM, "phyAM_27sp.tre")
# write.tree(phyEA, "phyEA_24sp.tre")


###################
# GLOBAL HL (48sp)
###################
hl <- cbind(morpho$sp, log(morpho$HL))
rownames(hl) <- hl[,1]
hl[,1] == tr55$tip.label
hl <- hl[tr55$tip.label,]
hl <- hl[,-1]
hl <- as.data.frame(hl)
hl <- na.omit(hl)
# write.table(hl, "hl.txt")

tr55
tr48 <- get_subtree_with_tips(tr55, rownames(hl))
tr48 <- tr48$subtree
# write.tree(tr48, "tr48_hl.tre")
rownames(hl) == tr48$tip.label  # CHECK

# Get priors
setBAMMpriors(phy = tr48, traits = "hl.txt", outfile = "myPriors_hl.txt")


##############
# EURASIA HL
##############
hlEA <- cbind(morphoEA$sp, log(morphoEA$HL))
rownames(hlEA) <- hlEA[,1]            
