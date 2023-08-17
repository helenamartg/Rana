##################
# Phylosignal
#################
rm(list=ls())
library(phytools)
library(ape)
library(geomorph)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

# 1. Import phylogeny
#=====================
phyloRana <- read.tree("DATA/phyloRana_74sp.tre")
phyloRana$tip.label

phylo68 <- read.tree("DATA/phylorana_68sp.tre")

# 2. Import traits
#===================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
head(dt74)
dt74 <- dt74[,-1]

dt68 <-  read.csv("DATA/dt_68sp.csv")
rownames(dt68) <-  dt68$sp
dt68 <- dt68[,-1]


# 3. Import climatic variables
#=============================
# clim <- read.csv("DATA/clim_74sp.csv")
# rownames(clim) <- clim$X # 74 observations 

clim_ea <- read.csv("DATA/clim_EA_31sp.csv")
rownames(clim_ea) <-  clim_ea$X
clim_am <- read.csv("DATA/clim_AM_37sp.csv")
rownames(clim_am) <- clim_am$X  # 68 observations

clim68 <- rbind(clim_ea, clim_am)


# 4. Check order
# #===============
# phyloRana$tip.label == dt74$sp
# dt74 <- dt74[phyloRana$tip.label,]

# phyloRana$tip.label == clim$X
# clim <- clim[phyloRana$tip.label,]

phylo68$tip.label == dt68$sp
dt68 <- dt68[phylo68$tip.label,]

phylo68$tip.label==clim68$X
clim68 <- clim68[phylo68$tip.label,]


# remove SD
clim2 <- clim68[,-grep("sd", colnames(clim68))]
clim2 <- clim2[,-c(1,2)]


# 5. phylogenetic signal  
#========================
# 5.1 SVL (68 sp)
#=================
svl <- setNames(dt68$log_svl, rownames(dt68))
physig_SVL_K <- phylosig(phylo68, svl, method = "K", test=TRUE)
physig_SVL_K
# Now we want to test for phylogenetic signal different than 1 (K=1 means species are as similar as expected given a Brownian motion model of evolution)

# First generate data with strong phylogenetic signal (null distribution for K=1)
nullK <- apply(fastBM(phylo68, n=1000, sig2=mean(pic(svl, phylo68)^2)),2,phylosig, tree = phylo68)
p <- mean(abs(log(c(physig_SVL_K$K, nullK)))>=abs(log(physig_SVL_K$K)))
round(p, 3)
# p = 0.066 which makes phylogenetic signal of SVL not significantly different from K = 1 
# There is phylogenetic signal


# 5.2 Climatic data (68 sp)
#===========================
clim2
physig_clim <- matrix(ncol = ncol(clim2), nrow = 5)
for ( i in 1:ncol(clim2)){
  sg <- setNames(clim2[,i], rownames(clim2))
  p <- phylosig(phylo68, sg, method = "lambda", test=TRUE)
  q <- phylosig(phylo68, sg, method = "K", test=TRUE)
  physig_clim[1,i] <- round(p$lambda, 3)
  physig_clim[2,i] <- p$P
  physig_clim[3,i] <- round(q$K, 3)
  physig_clim[4,i] <- q$P
  nullK <- apply(fastBM(phylo68, n=1000, sig2=mean(pic(sg, phylo68)^2)),2, phylosig, tree = phylo68)
  p <- mean(abs(log(c(q$K, nullK)))>=abs(log(q$K)))
  physig_clim[5,i] <- round(p, 3)
  
}

rownames(physig_clim) <- c("lambda", "p-valor", "Blombergs'K", "p-valor", "K not 1 (p-value)")
colnames(physig_clim) <- colnames(clim2)
x <- t(physig_clim)

#SAVE
# write.table(x, "physignal_clim.csv")



#=============================
# 6. PHYLOSIGNAL by radiation
#=============================
dt68 <- read.csv("DATA/dt_68sp.csv")
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

svlEA <- setNames(dtEA$log_svl, dtEA$sp)
svlAM <- setNames(dtAM$log_svl, dtAM$sp)

tr_EA <- read.tree("DATA/phyloEA_31sp.tre")
tr_AM <- read.tree("DATA/phyloAM_37sp.tre")


# 6.1 SVL
#=========
# EURASIA (31sp)
physig_SVL_EA <- phylosig(tr_EA, svlEA, method = "K", test=TRUE)  # Significant
nullK <- apply(fastBM(tr_EA, n=1000, sig2=mean(pic(svlEA, tr_EA)^2)),2,phylosig, tree = tr_EA)
p <- mean(abs(log(c(physig_SVL_EA$K, nullK)))>=abs(log(physig_SVL_EA$K)))
round(p, 3)  # Not significantly different from 1 (phylogenetic signal)

# AMERICA (37sp)
physig_SVL_AM <- phylosig(tr_AM, svlAM, method = "K", test=TRUE)  # Not significant
nullK <- apply(fastBM(tr_AM, n=1000, sig2=mean(pic(svlAM, tr_AM)^2)),2,phylosig, tree = tr_AM)
p <- mean(abs(log(c(physig_SVL_AM$K, nullK)))>=abs(log(physig_SVL_AM$K)))
round(p, 3)  #Significantly different from 1 (NO phyogenetic signal)


# 6.2 CLIMATE 
#==============
# EURASIA
clim_EA <- read.csv("DATA/clim_EA_31sp.csv")
rownames(clim_EA) <- clim_EA$X
clim_EA <- clim_EA[,-c(1,2)]
clim_EA <- clim_EA[,-grep("sd", colnames(clim_EA))]

physig_clim <- matrix(ncol = ncol(clim_EA), nrow = 5)
for ( i in 1:ncol(clim_EA)){
  sg <- setNames(clim_EA[,i], rownames(clim_EA))
  p <- phylosig(tr_EA, sg, method = "lambda", test=TRUE)
  q <- phylosig(tr_EA, sg, method = "K", test=TRUE)
  physig_clim[1,i] <- round(p$lambda, 3)
  physig_clim[2,i] <- p$P
  physig_clim[3,i] <- round(q$K, 3)
  physig_clim[4,i] <- q$P
  nullK <- apply(fastBM(tr_EA, n=1000, sig2=mean(pic(sg, tr_EA)^2)),2, phylosig, tree = tr_EA)
  p <- mean(abs(log(c(q$K, nullK)))>=abs(log(q$K)))
  physig_clim[5,i] <- round(p, 3)
  
}

rownames(physig_clim) <- c("lambda", "p-valor", "Blombergs'K", "p-valor", "K not 1 (p-value)")
colnames(physig_clim) <- colnames(clim2)
x <- t(physig_clim)

#SAVE
# write.csv(x, "phylo_clim_EA.csv")

# AMERICA
clim_AM <- read.csv("DATA/clim_AM_37sp.csv")
rownames(clim_AM) <- clim_AM$X
clim_AM <- clim_AM[,-c(1,2)]
clim_AM <- clim_AM[,-grep("sd", colnames(clim_AM))]

physig_clim <- matrix(ncol = ncol(clim_AM), nrow = 5)
for ( i in 1:ncol(clim_AM)){
  sg <- setNames(clim_AM[,i], rownames(clim_AM))
  p <- phylosig(tr_AM, sg, method = "lambda", test=TRUE)
  q <- phylosig(tr_AM, sg, method = "K", test=TRUE)
  physig_clim[1,i] <- round(p$lambda, 3)
  physig_clim[2,i] <- p$P
  physig_clim[3,i] <- round(q$K, 3)
  physig_clim[4,i] <- q$P
  nullK <- apply(fastBM(tr_AM, n=1000, sig2=mean(pic(sg, tr_AM)^2)),2, phylosig, tree = tr_AM)
  p <- mean(abs(log(c(q$K, nullK)))>=abs(log(q$K)))
  physig_clim[5,i] <- round(p, 3)
  
}

rownames(physig_clim) <- c("lambda", "p-valor", "Blombergs'K", "p-valor", "K not 1 (p-value)")
colnames(physig_clim) <- colnames(clim2)
x <- t(physig_clim)

#SAVE
# write.csv(x, "phylo_clim_AM.csv")


# 7. Phylogenetic signal of morphological variables (GLOBAL and by radiation)
#==================================================
library(castor)
morpho <- read.csv("DATA/dt_morpho_55sp.csv")
rownames(morpho) <- morpho$sp

# Extract 6 species from EA radiation that are distributed in aM
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]

# 4 morphological variables: HL, HW, TL, FL with lack of data for some species
# We calculate phylogenetic signal for each variable, creating a loop to clip phylogeny in each case
# because sample size is different for each variable.

# Extract columns with residuals only
morpho <- morpho[,c(1,4,14:17)]
morpho <- na.omit(morpho)  # 34sp

morpho_EA <- morpho[morpho$phylo_reg=="EA",] # 21sp
morpho_AM <- morpho[morpho$phylo_reg=="AM",] # 13 sp

# clip phylogeny
phylo34 <- get_subtree_with_tips(phylo68, morpho$X)
phylo34 <- phylo34$subtree

phylo_EA <- get_subtree_with_tips(phylo34, morpho_EA$X)  # Remove sp from phylogeny
phylo_EA <- phylo_EA$subtree 

phylo_AM <- get_subtree_with_tips(phylo34, morpho_AM$X)  # Remove sp from phylogeny
phylo_AM <- phylo_AM$subtree 


# LOOP for all (34 sp)
#==============
mat <- NULL
for (i in 3:ncol(morpho)){
  # General
  x <- setNames(morpho[,i], morpho$X)   
  
  physig_x <- phylosig(phylo34, x, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo34, n=1000, 
                        sig2=mean(pic(x, phylo34)^2)),2,phylosig, tree = phylo34) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x $K, nullK)))>=abs(log(physig_x $K)))
  round(p, 3)  
  
  temp <- c(physig_x$K, physig_x$P, p)
  mat <- rbind(mat, temp)
  
  # Eurasia
  x_EA <- morpho[morpho$phylo_reg=="EA",]  # Separate Eurasian radiation
  
  x_EA <- setNames(x_EA[,i], x_EA$X)
  
  physig_x_EA <- phylosig(phylo_EA, x_EA, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo_EA, n=1000, 
                        sig2=mean(pic(x_EA, phylo_EA)^2)),2,phylosig, tree = phylo_EA) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x_EA$K, nullK)))>=abs(log(physig_x_EA$K)))
  round(p, 3) 
  
  temp.ea <- c(physig_x_EA$K, physig_x_EA$P, p)
  mat <- rbind(mat, temp.ea)
  
  
  # America
  x_AM <- morpho[morpho$phylo_reg=="AM",] # Separate American radiation
  
  x_AM <- setNames(x_AM[,i], x_AM$X)
  
  physig_x_AM <- phylosig(phylo_AM, x_AM, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo_AM, n=1000, 
                        sig2=mean(pic(x_AM, phylo_AM)^2)),2,phylosig, tree = phylo_AM) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x_AM$K, nullK)))>=abs(log(physig_x_AM$K)))
  round(p, 3) 
  
  temp.am <- c(physig_x_AM$K, physig_x_AM$P, p)
  mat <- rbind(mat, temp.am)
  
}


colnames(mat) <- c("Blombergs'K", "P-value", "K not 1 (pvalue)")
rownames(mat) <-  c("Global_HL", "Eurasian_HL", "American_HL", "Global_HW", "Eurasian_HW", "American_HW",
                    "GLobal_FL", "Eurasian_FL", "American_FL", "Global_TL", "Eurasian_TL", "American_TL")

# SAVE
# write.csv(mat, "phylosig_morpho.csv")




###########################################
# clipping phylogeny
##########################################

# LOOP for all (clipping phylogeny)
#==============
mat <- NULL
for (i in 3:ncol(morpho)){
  # General
  morpho.x <- morpho[complete.cases(morpho[i]),]  # Remove sp with NAs
  phylo.x <- get_subtree_with_tips(phyloRana, morpho.x$sp)  # Remove sp from phylogeny
  phylo.x <- phylo.x$subtree
  
  x <- setNames(morpho.x[,i], morpho.x$sp)   
  
  physig_x <- phylosig(phylo.x, x, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.x, n=1000, 
                        sig2=mean(pic(x, phylo.x)^2)),2,phylosig, tree = phylo.x) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x $K, nullK)))>=abs(log(physig_x $K)))
  round(p, 3)  
  
  temp <- c(physig_x$K, physig_x$P, p)
  mat <- rbind(mat, temp)
  
  # Eurasia
  x_EA <- morpho.x[morpho.x$phylo_reg=="EA",]  # Separate Eurasian radiation
  
  phylo.x.EA <- get_subtree_with_tips(phylo.x, x_EA$sp)  # Remove sp from phylogeny
  phylo.x.EA <- phylo.x.EA$subtree  
  
  x_EA <- setNames(x_EA[,i], x_EA$sp)
  
  physig_x_EA <- phylosig(phylo.x.EA, x_EA, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.x.EA, n=1000, 
                        sig2=mean(pic(x_EA, phylo.x.EA)^2)),2,phylosig, tree = phylo.x.EA) # Null distribution of K=1
  p <- mean(abs(log(c(physig_x_EA$K, nullK)))>=abs(log(physig_x_EA$K)))
  round(p, 3) 
  
  temp.ea <- c(physig_x_EA$K, physig_x_EA$P, p)
  mat <- rbind(mat, temp.ea)
  
  
  # America
  x_AM <- morpho.x[morpho.x$phylo_reg=="AM",] # Separate American radiation
  
  phylo.x.AM <- get_subtree_with_tips(phylo.x, x_AM$sp) # Remove sp from phylogeny
  phylo.x.AM <- phylo.x.AM$subtree 
  
  x_AM <- setNames(x_AM[,i], x_AM$sp)
  
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
# write.csv(mat, "phylosig_morpho_loop.csv")

