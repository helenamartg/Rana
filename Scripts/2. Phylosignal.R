##################
# Phylosignal
#################
rm(list=ls())
library(phytools)
library(ape)
library(geomorph)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4publication/2nd_Revision/Reanalisis")

# 1. Import phylogeny
#=====================
phylo59 <- read.tree("DATA/phyloRana_59sp.tre")
phylo59$tip.label

phylo45 <- read.tree("DATA/phylorana_45sp.tre")


# 2. Import traits
#===================
dt59 <- read.csv("DATA/dt_59sp.csv", header=T)
rownames(dt59) <- dt59$Species
head(dt59)
dt59 <- dt59[,-1]


# 4. Check order
# #===============
phylo59$tip.label == dt59$Species
dt59 <- dt59[phylo59$tip.label,]



# 5. phylogenetic signal  
#========================
# 5.1 SVL (59 sp)
#=================
svl <- setNames(dt59$log_svl, rownames(dt59))
physig_SVL_K <- phylosig(phylo59, svl, method = "K", test=TRUE)
physig_SVL_K
# Now we want to test for phylogenetic signal different than 1 (K=1 means species are as similar as expected given a Brownian motion model of evolution)

# First generate data with strong phylogenetic signal (null distribution for K=1)
nullK <- apply(fastBM(phylo59, n=1000, sig2=mean(pic(svl, phylo59)^2)),2,phylosig, tree = phylo59)
p <- mean(abs(log(c(physig_SVL_K$K, nullK)))>=abs(log(physig_SVL_K$K)))
round(p, 3)
# p = 0.069 which makes phylogenetic signal of SVL not significantly different from K = 1 
# There is phylogenetic signal



#=============================
# 6. PHYLOSIGNAL by radiation
#=============================
dtEA <- dt59[dt59$real_reg=="EA",]
dtAM <- dt59[dt59$real_reg=="AM",]

tr_EA <- read.tree("Data/phyloEA_26sp.tre")
tr_AM <- read.tree("Data/phyloAM_33sp.tre")

# Check order
dtEA$Species==tr_EA$tip.label
dtAM$Species==tr_AM$tip.label

svlEA <- setNames(dtEA$log_svl, dtEA$Species)
svlAM <- setNames(dtAM$log_svl, dtAM$Species)


# 6.1 SVL
#=========
# EURASIA (26sp)
physig_SVL_EA <- phylosig(tr_EA, svlEA, method = "K", test=TRUE)  # Significant
nullK <- apply(fastBM(tr_EA, n=1000, sig2=mean(pic(svlEA, tr_EA)^2)),2,phylosig, tree = tr_EA)
p <- mean(abs(log(c(physig_SVL_EA$K, nullK)))>=abs(log(physig_SVL_EA$K)))
round(p, 3)  # Not significantly different from 1 (phylogenetic signal)

# AMERICA (33sp)
physig_SVL_AM <- phylosig(tr_AM, svlAM, method = "K", test=TRUE)  # Not significant
nullK <- apply(fastBM(tr_AM, n=1000, sig2=mean(pic(svlAM, tr_AM)^2)),2,phylosig, tree = tr_AM)
p <- mean(abs(log(c(physig_SVL_AM$K, nullK)))>=abs(log(physig_SVL_AM$K)))
round(p, 3)  #Significantly different from 1 (NO phyogenetic signal)



# 7. Phylogenetic signal of morphological variables (GLOBAL and by radiation)
#==================================================
library(castor)
morpho <- read.csv("Data/morpho_phylores_45sp.csv")
rownames(morpho) <- morpho$Species

# Extract 6 species from EA radiation that are distributed in aM
morpho <- morpho[morpho$Species!="Rana_tlaloci",]
morpho <- morpho[morpho$Species!="Rana_luteiventris",]
morpho <- morpho[morpho$Species!="Rana_boylii",]
morpho <- morpho[morpho$Species!="Rana_sierrae",]
morpho <- morpho[morpho$Species!="Rana_muscosa",]
morpho <- morpho[morpho$Species!="Rana_aurora",]
morpho <- morpho[morpho$Species!="Rana_cascadae",]

# 4 morphological variables: HL, HW, TL, FL with lack of data for some species
# We calculate phylogenetic signal for each variable, creating a loop to clip phylogeny in each case
# because sample size is different for each variable.

# Extract columns with residuals only
morpho <- morpho[,c(2,4,13:16)]
na.omit(morpho)  # 26sp

# morpho_EA <- morpho[morpho$phylo_reg=="EA",] # 16sp
# morpho_AM <- morpho[morpho$phylo_reg=="AM",] # 10 sp
# 
# # clip phylogeny
# phylo <- get_subtree_with_tips(phylo59, morpho$X)
# phylo34 <- phylo34$subtree
# 
# phylo_EA <- get_subtree_with_tips(phylo34, morpho_EA$X)  # Remove sp from phylogeny
# phylo_EA <- phylo_EA$subtree 
# 
# phylo_AM <- get_subtree_with_tips(phylo34, morpho_AM$X)  # Remove sp from phylogeny
# phylo_AM <- phylo_AM$subtree 



# LOOP for all (clipping phylogeny)
#==============
mat <- NULL
for (i in 3:ncol(morpho)){
  # General
  morpho.x <- morpho[complete.cases(morpho[i]),]  # Remove sp with NAs
  phylo.x <- get_subtree_with_tips(phylo59, morpho.x$Species)  # Remove sp from phylogeny
  phylo.x <- phylo.x$subtree
  
  x <- setNames(morpho.x[,i], morpho.x$Species)   
  
  physig_x <- phylosig(phylo.x, x, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.x, n=1000, 
                        sig2=mean(pic(x, phylo.x)^2)),2,phylosig, tree = phylo.x) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x $K, nullK)))>=abs(log(physig_x $K)))
  round(p, 3)  
  
  temp <- c(physig_x$K, physig_x$P, p)
  mat <- rbind(mat, temp)
  
  # Eurasia
  x_EA <- morpho.x[morpho.x$phylo_reg=="EA",]  # Separate Eurasian radiation
  
  phylo.x.EA <- get_subtree_with_tips(phylo.x, x_EA$Species)  # Remove sp from phylogeny
  phylo.x.EA <- phylo.x.EA$subtree  
  
  x_EA <- setNames(x_EA[,i], x_EA$Species)
  
  physig_x_EA <- phylosig(phylo.x.EA, x_EA, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.x.EA, n=1000, 
                        sig2=mean(pic(x_EA, phylo.x.EA)^2)),2,phylosig, tree = phylo.x.EA) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x_EA$K, nullK)))>=abs(log(physig_x_EA$K)))
  round(p, 3) 
  
  temp.ea <- c(physig_x_EA$K, physig_x_EA$P, p)
  mat <- rbind(mat, temp.ea)
  
  
  # America
  x_AM <- morpho.x[morpho.x$phylo_reg=="AM",] # Separate American radiation
  
  phylo.x.AM <- get_subtree_with_tips(phylo.x, x_AM$Species) # Remove sp from phylogeny
  phylo.x.AM <- phylo.x.AM$subtree 
  
  x_AM <- setNames(x_AM[,i], x_AM$Species)
  
  physig_x_AM <- phylosig(phylo.x.AM, x_AM, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.x.AM, n=1000, 
                        sig2=mean(pic(x_AM, phylo.x.AM)^2)),2,phylosig, tree = phylo.x.AM) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x_AM$K, nullK)))>=abs(log(physig_x_AM$K)))
  round(p, 3) 
  
  temp.am <- c(physig_x_AM$K, physig_x_AM$P, p)
  mat <- rbind(mat, temp.am)
  
}


colnames(mat) <- c("Blombergs'K", "P-value", "K not 1 (pvalue)")
rownames(mat) <-  c("Global_HL", "Eurasian_HL", "American_HL", "Global_HW", "Eurasian_HW", "American_HW",
                    "GLobal_FL", "Eurasian_FL", "American_FL", "Global_TL", "Eurasian_TL", "American_TL")

# SAVE
# write.csv(mat, "Results/phylosig_morpho_loop.csv")

