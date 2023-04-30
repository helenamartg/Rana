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

# 2. Import traits
#===================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
head(dt74)
dt74 <- dt74[,-1]


# 3. Import climatic variables
#=============================

clim <- read.csv("DATA/clim_74sp.csv")
rownames(clim) <- clim$X
# 74 observations 


# 4. Check order
#===============
phyloRana$tip.label == dt74$sp
dt74 <- dt74[phyloRana$tip.label,]

phyloRana$tip.label == clim$X
clim <- clim[phyloRana$tip.label,]

# remove SD
clim2 <- clim[,-grep("sd", colnames(clim))]
clim2 <- clim2[,-c(1,2)]


# 5. phylogenetic signal  
#========================

# 5.1 SVL (74 sp)
#=================
svl <- setNames(dt74$SVL, rownames(dt74))
physig_SVL_K <- phylosig(phyloRana, svl, method = "K", test=TRUE)
physig_SVL_K
physig_SVL_lambda <- phylosig(phyloRana, svl, method = "lambda", test = T) #lambda si sale significativo
physig_SVL_lambda

# Now we want to test for phylogenetic signal different than 1 (K=1 means species are as similar as expected given a Brownian motion model of evolution)

# First generate data with strong phylogenetic signal (null distribution for K=1)
nullK <- apply(fastBM(phyloRana, n=1000, sig2=mean(pic(svl, phyloRana)^2)),2,phylosig, tree = phyloRana)
p <- mean(abs(log(c(physig_SVL_K$K, nullK)))>=abs(log(physig_SVL_K$K)))
round(p, 3)
# p = 0.003 which makes phylogenetic signal of SVL significantly different from K = 1 
# There is no phylogenetic signal


# 5.2 Climatic data (74 sp)
#===========================
clim2
physig_clim <- matrix(ncol = ncol(clim2), nrow = 5)
for ( i in 1:ncol(clim2)){
  sg <- setNames(clim2[,i], rownames(clim2))
  p <- phylosig(phyloRana, sg, method = "lambda", test=TRUE)
  q <- phylosig(phyloRana, sg, method = "K", test=TRUE)
  physig_clim[1,i] <- round(p$lambda, 3)
  physig_clim[2,i] <- p$P
  physig_clim[3,i] <- round(q$K, 3)
  physig_clim[4,i] <- q$P
  nullK <- apply(fastBM(phyloRana, n=1000, sig2=mean(pic(sg, phyloRana)^2)),2, phylosig, tree = phyloRana)
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
dt68 <- dt68[,-1]
rownames(dt68) <- dt68$sp
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

svlEA <- setNames(dtEA$SVL, dtEA$sp)
svlAM <- setNames(dtAM$SVL, dtAM$sp)

tr_EA <- read.tree("DATA/phyloEA_31sp.tre")
tr_AM <- read.tree("DATA/phyloAM_37sp.tre")


# 6.1 SVL
#=========
# EURASIA
physig_SVL_EA <- phylosig(tr_EA, svlEA, method = "K", test=TRUE)  # Significant
nullK <- apply(fastBM(tr_EA, n=1000, sig2=mean(pic(svlEA, tr_EA)^2)),2,phylosig, tree = tr_EA)
p <- mean(abs(log(c(physig_SVL_EA$K, nullK)))>=abs(log(physig_SVL_EA$K)))
round(p, 3)  # Not significantly different from 1 (phylogenetic signal)

# AMERICA
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

# Extract 6 species from EA radiation that are distributed in aM
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]

#4 morphological variables: HL, HW, TL, FL with lack of data for some species
# We calculate phylogenetic signal for each variable, creating a loop to clip phylogeny in each case
# because sample size is different for each variable.

# LOOP for all
#==============
mat <- NULL
for (i in 9:ncol(morpho)){
  # General
  morpho.x <- morpho[complete.cases(morpho[i]),]  # Remove sp with NAs
  phylo.x <- get_subtree_with_tips(phyloRana, morpho.x$sp)  # Remove sp from phylogeny
  phylo.x <- phylo.x$subtree
  
  res <- resid(lm(morpho.x[,i]~morpho.x$SVL)) # Residuals
  
  x <- setNames(res, morpho.x$sp)   
  
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
  
  res.ea <- resid(lm(x_EA[,i]~x_EA$SVL)) # Residuals
  
  x_EA <- setNames(res.ea, x_EA$sp)
  
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
  
  res.am <- resid(lm(x_AM[,i]~x_AM$SVL)) # Residuals
  
  x_AM <- setNames(res.am, x_AM$sp)
  
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
# write.csv(mat, "phylogenetic_signal_morpho_RESIDUALS.csv")



###################
# REDUCED DATASET 
####################
# Reduce the dataset to species that have values for all morphological variables (same sample size)
morpho.red <- na.omit(morpho[,c(2:5,8:12)])
table(morpho.red$phylo_reg)

# Calculate Residuals
mat.red <- NULL
for (i in 6:ncol(morpho.red)){
  res <- resid(lm(log(morpho.red[,i]) ~ log(morpho.red$SVL)))
  mat.red <- cbind(mat.red, res)
}

colnames(mat.red) <- c("res_HL", "res_HW", "res_FL", "res_TL")
rownames(mat.red) <- morpho.red$sp

# Merge datasets
mat.red <- cbind(morpho.red, mat.red) # 34 sp

# Clip phylogeny
phylo.red <- get_subtree_with_tips(phylo.x, mat.red$sp)
phylo.red <- phylo.red$subtree  # 34 sp

# Separate by radiations
mat.red.EA <- mat.red[mat.red$phylo_reg=="EA",]
mat.red.AM <- mat.red[mat.red$phylo_reg=="AM",]

phylo.red.EA <- (get_subtree_with_tips(phylo.x, mat.red.EA$sp))$subtree  # 21 sp
phylo.red.AM <- (get_subtree_with_tips(phylo.x, mat.red.AM$sp))$subtree  # 13 sp


# PHYLOG. SIGNAL GLOBAL
temp <- NULL
m <- NULL
for (i in 10:ncol(mat.red)){
  x <- setNames(mat.red[,i], mat.red$sp)
  
  physig <- phylosig(phylo.red, x, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.red, n=1000, 
                        sig2=mean(pic(x, phylo.red)^2)),2,phylosig, tree = phylo.red) # Null distribution of K=1
  p <- mean(abs(log(c(physig$K, nullK)))>=abs(log(physig$K)))
  round(p, 3) 
  
  temp <- c(physig$K, physig$P, p)
  m <- rbind(m, temp)
}

colnames(m) <- c("Blombergs'K", "P-value", "K not 1 (pvalue)")
rownames(m) <- c("Global_res_HL", "Global_res_HW", "Global_res_FL", "Global_res_TL")



# PHYLOSIGNAL EURASIA
temp.ea <- NULL
m.ea <- NULL
for (i in 10:ncol(mat.red.EA)){
  x <- setNames(mat.red.EA[,i], mat.red.EA$sp)
  
  physig <- phylosig(phylo.red.EA, x, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.red.EA, n=1000, 
                        sig2=mean(pic(x, phylo.red.EA)^2)),2,phylosig, tree = phylo.red.EA) # Null distribution of K=1
  p <- mean(abs(log(c(physig$K, nullK)))>=abs(log(physig$K)))
  round(p, 3) 
  
  temp.ea <- c(physig$K, physig$P, p)
  m.ea <- rbind(m.ea, temp.ea)
}

colnames(m.ea) <- c("Blombergs'K", "P-value", "K not 1 (pvalue)")
rownames(m.ea) <- c("Eurasia_res_HL", "Eurasia_res_HW", "Eurasia_res_FL", "Eurasia_res_TL")


# PHYLOSIGNAL AMERICA
temp.am <- NULL
m.am <- NULL
for (i in 10:ncol(mat.red.AM)){
  x <- setNames(mat.red.AM[,i], mat.red.AM$sp)
  
  physig <- phylosig(phylo.red.AM, x, method = "K", test=TRUE) 
  nullK <- apply(fastBM(phylo.red.AM, n=1000, 
                        sig2=mean(pic(x, phylo.red.AM)^2)),2,phylosig, tree = phylo.red.AM) # Null distribution of K=1
  p <- mean(abs(log(c(physig$K, nullK)))>=abs(log(physig$K)))
  round(p, 3) 
  
  temp.am <- c(physig$K, physig$P, p)
  m.am <- rbind(m.am, temp.am)
}

colnames(m.am) <- c("Blombergs'K", "P-value", "K not 1 (pvalue)")
rownames(m.am) <- c("America_res_HL", "America_res_HW", "America_res_FL", "America_res_TL")


# Merge the three datasets
merged <- rbind(m, m.ea, m.am)

# SAVE
# write.csv(merged, "physig_morpho_reduced.csv")
