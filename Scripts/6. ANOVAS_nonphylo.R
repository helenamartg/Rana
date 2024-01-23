##################################
#  NEW ANOVAS (non phylogenetic)
##################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(RRPP)
library(ape)
library(phytools)
library(castor)

# 1. Import data
#================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
dt74 <- dt74[,-1]

dt68 <-read.csv("DATA/dt_68sp.csv")
dt68 <- dt68[,-1]
rownames(dt68) <- dt68$sp

# 2. Import phylogeny
#=====================
phyloRana <- read.tree("DATA/phyloRana_74sp.tre")
phyloRana$tip.label
phyloRana_68 <- read.tree("DATA/phyloRana_68sp.tre")

# 3. Match with phylogeny
#=========================
dt74 <- dt74[phyloRana$tip.label,]
dt68 <-  dt68[phyloRana_68$tip.label,]

dt74$sp == phyloRana$tip.label
dt68$sp == phyloRana_68$tip.label

# 4. ANOVAs
#===========
rrpp.global <- rrpp.data.frame(svl = dt74$log_svl, 
                               sp = as.factor(dt74$sp),
                               radiation = as.factor(dt74$phylo_reg))

rrpp.global68 <- rrpp.data.frame(svl = dt68$log_svl, 
                                 sp = as.factor(dt68$sp),
                                 radiation = as.factor(dt68$phylo_reg))

# Are body sizes different across radiations?
lm.SVL <- lm.rrpp(svl ~ radiation, data = rrpp.global68, SS.type = "III")  # Without phylogeny (GLS)
anova(lm.SVL)   # Significant


# MORPHO
#========
# Remaining morphological variables
morpho <- read.csv("DATA/dt_morpho_55sp.csv")
rownames(morpho) <- morpho$sp

#remove 6 species from Eurasia
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!= "Rana_tlaloci",]

morpho.res <- morpho[,c(2,4,14:17)]

# MEGALOOP to not miss species and clip dataset and phylogeny by each trait
mat <- NULL
for (i in 3:6){
  morpho.x <- morpho.res[complete.cases(morpho.res[,i]),]      # Remove sp from morpho dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.x$sp))$subtree  # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  
  cov.temp <- vcv(phylo.x, model="Brownian")
  
  rrpp.morpho <- rrpp.data.frame(morpho = morpho.x[,i],
                                 sp = as.factor(morpho.x$sp),
                                 radiation = as.factor(morpho.x$phylo_reg))
  
  lm.rad <- lm.rrpp(morpho ~ radiation, data = rrpp.morpho, SS.type = "III")
  anova(lm.rad)
  
  #MERGE
  names(morpho.x)
  lms <- cbind(anova(lm.rad)$table, names(morpho.x[i]))
  
  mat <- rbind(mat, lms)
  
}


# Merge SVL and morpho
svl <- cbind(anova(lm.SVL)$table, "log_svl")
names(svl) <- names(mat)
mat

def <- rbind(svl, mat)


# SAVE
# write.table(def, "results/Anovas_nonphylogenetic.csv", sep=";")



# PLOT
layout(matrix(c(1,2,3,4),nrow=2, ncol=2, byrow = T))
par(mar=c(4,4,4,4))
boxplot(morpho$res_HL ~ morpho$phylo_reg, ylab="res_HL", xlab="")
boxplot(morpho$res_HW ~ morpho$phylo_reg, ylab="res_HW", xlab="")
boxplot(morpho$res_FL ~ morpho$phylo_reg, ylab="res_FL", xlab="")
boxplot(morpho$res_TL ~ morpho$phylo_reg, ylab="res_TL", xlab="")

