#################################
# Phylogenetic linear models
# with BM and OU: 
# MORPHOLOGICAL VARIABLES
#################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(phylolm)
library(phytools)


# 1. Import traits: morphological traits 
#===========================================
morpho <- read.csv("DATA/dt_phylores_55sp.csv")
morpho <- morpho[,-1]
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]
rownames(morpho) <- morpho$sp

# 2. Import phylogeny
#======================
phyloEA <- read.tree("DATA/phyloEA_31sp.tre")
phyloAM <- read.tree("DATA/phyloAM_37sp.tre")


# 3. Import climatic variables
#===============================
climEA <- read.csv("DATA/Clim_EA_31sp.csv")
climAM <- read.csv("DATA/Clim_AM_37sp.csv")

climEA <- climEA[,-1]
climAM <- climAM[,-1]
rownames(climEA) <- climEA$X
rownames(climAM) <- climAM$X

# Scale climatic variables
climEA <- climEA[,-grep("SD", colnames(climEA))]
climEA <- scale(climEA[,-1])

climAM <- climAM[,-grep("SD", colnames(climAM))]
climAM <- scale(climAM[,-1])


# 4. Check order
#================
rownames(climEA)==phyloEA$tip.label
rownames(climAM)==phyloAM$tip.label


# 6.1. EURASIA
#=================
morpho.ea <- morpho[morpho$phylo_reg=="EA",]
morpho.ea <- morpho.ea[, c(1,13:16)]


results_morpho_ea <- data.frame()
mat.bm.ea <- data.frame()
mat.ou.ea <- data.frame()

for(i in 2:ncol(morpho.ea)) {
  morpho.x <- morpho.ea[complete.cases(morpho.ea[,i]), c(1,i)]          # Remove sp from morpho dataset
  clim.x <- climEA[morpho.ea$sp,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloEA, morpho.x$sp))$subtree      # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- clim.x[phylo.x$tip.label,]
  
  trait <- setNames(morpho.x[,2], morpho.x$sp)
  
  for (j in 1:ncol(clim.x)) {
    
    bm_model_eurasia <- phylolm(trait ~ clim.x[,j], phy = phylo.x, model = "BM")
    ou_model_eurasia <- phylolm(trait ~ clim.x[,j], phy = phylo.x, model = "OUfixedRoot")
    
    bm_summary <- summary(bm_model_eurasia)
    ou_summary <- summary(ou_model_eurasia)
    
    bm <- as.data.frame(cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik))
    bm$trait <- as.character(colnames(morpho.x)[2])
    rownames(bm)[2] <- colnames(clim.x)[j]
    
    ou <- as.data.frame(cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik))
    ou$trait <- as.character(colnames(morpho.x)[2])
    rownames(ou)[2] <- colnames(clim.x)[j]
    
    mat.bm.ea <- rbind(mat.bm.ea, bm)
    mat.ou.ea <- rbind(mat.ou.ea, ou)
    
    AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
    AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)
    
    x <- cbind(trait = colnames(morpho.x)[2], 
               clim_var = colnames(clim.x)[j], AICc_BM = as.numeric(AICc_bm), AICc_OU = as.numeric(AICc_ou),
               Radiation = "Eurasia", Lowest_AICc=NA)
    
    results_morpho_ea <- rbind(results_morpho_ea, x)
    
  }
}

results_morpho_ea <- as.data.frame(results_morpho_ea)
results_morpho_ea[,3] <- as.numeric(results_morpho_ea[,3])
results_morpho_ea[,4] <- as.numeric(results_morpho_ea[,4])

# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results_morpho_ea)){
  if (results_morpho_ea[i,]$AICc_BM < results_morpho_ea[i,]$AICc_OU) {
    results_morpho_ea[i,]$Lowest_AICc <- "BM model"
  } else {
    results_morpho_ea[i,]$Lowest_AICc <- "OU model"
  }
}


# Calculate delta AICc
results_morpho_ea <- results_morpho_ea %>%
  mutate(
    Delta_AIC_BM = AICc_BM - pmin(AICc_BM, AICc_OU),
    Delta_AIC_OU = AICc_OU - pmin(AICc_BM, AICc_OU)
  )

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 and 4 units
results_morpho_ea <- results_morpho_ea %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

results_morpho_ea <- results_morpho_ea %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)


# Select the best model 
results_morpho_ea$Best_model <- NA
for (i in 1:nrow(results_morpho_ea)){
  if (results_morpho_ea[i,10]==TRUE){
    results_morpho_ea[i,11] <- results_morpho_ea[i,6]
  } else {
    results_morpho_ea[i,11] <- "BM"
  }
}

print(results_morpho_ea)



# 6.2. AMERICA
#=================
library(ape)
morpho.am <- morpho[morpho$phylo_reg=="AM",]
morpho.am <- morpho.am[, c(1,13:16)]


results_morpho_am <- data.frame()
mat.bm.am <- data.frame()
mat.ou.am <- data.frame()

for(i in 2:ncol(morpho.am)) {
  morpho.x <- morpho.am[complete.cases(morpho.am[,i]), c(1,i)]          # Remove sp from morpho dataset
  clim.x <- climAM[morpho.am$sp,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloAM, morpho.x$sp))$subtree      # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- clim.x[phylo.x$tip.label,]
  
  trait <- setNames(morpho.x[,2], morpho.x$sp)
  
  
  # Set the root edge to 0 (or remove it entirely)
  if(!is.null(phylo.x$root.edge)) {
    phylo.x$root.edge <- NULL  # Remove root edge if it exists
  }
  
  for (j in 1:ncol(clim.x)) {
    
    bm_model_america <- phylolm(trait ~ clim.x[,j], phy = phylo.x, model = "BM")
    ou_model_america <- phylolm(trait ~ clim.x[,j], phy = phylo.x, model = "OUfixedRoot")
    
    bm_summary <- summary(bm_model_america)
    ou_summary <- summary(ou_model_america)
    
    bm <- as.data.frame(cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik))
    bm$trait <- as.character(colnames(morpho.x)[2])
    rownames(bm)[2] <- colnames(clim.x)[j]
    
    ou <- as.data.frame(cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik))
    ou$trait <- as.character(colnames(morpho.x)[2])
    rownames(ou)[2] <- colnames(clim.x)[j]
    
    mat.bm.am <- rbind(mat.bm.am, bm)
    mat.ou.am <- rbind(mat.ou.am, ou)
    
    AICc_bm <- AIC(bm_model_america, corrected = TRUE)
    AICc_ou <- AIC(ou_model_america, corrected = TRUE)
    
    x <- cbind(trait = colnames(morpho.x)[2], 
               clim_var = colnames(clim.x)[j], AICc_BM = as.numeric(AICc_bm), AICc_OU = as.numeric(AICc_ou),
               Radiation = "America", Lowest_AICc=NA)
    
    results_morpho_am <- rbind(results_morpho_am, x)
    
  }
}


results_morpho_am <- as.data.frame(results_morpho_am)
results_morpho_am[,3] <- as.numeric(results_morpho_am[,3])
results_morpho_am[,4] <- as.numeric(results_morpho_am[,4])


# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results_morpho_am)){
  if (results_morpho_am[i,]$AICc_BM < results_morpho_am[i,]$AICc_OU) {
    results_morpho_am[i,]$Lowest_AICc <- "BM model"
  } else {
    results_morpho_am[i,]$Lowest_AICc <- "OU model"
  }
}


# Calculate delta AICc
results_morpho_am <- results_morpho_am %>%
  mutate(
    Delta_AIC_BM = AICc_BM - pmin(AICc_BM, AICc_OU),
    Delta_AIC_OU = AICc_OU - pmin(AICc_BM, AICc_OU)
  )

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 and 4 units
results_morpho_am <- results_morpho_am %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

results_morpho_am <- results_morpho_am %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)


# Select the best model 
results_morpho_am$Best_model <- NA
for (i in 1:nrow(results_morpho_am)){
  if (results_morpho_am[i,10]==TRUE){
    results_morpho_am[i,11] <- results_morpho_am[i,6]
  } else {
    results_morpho_am[i,11] <- "BM"
  }
}

print(results_morpho_am)



# SAVE RESULTS
################
# write.csv(results_morpho_ea, "Results/phylolm/phylolm_morpho_EA.csv")
# write.csv(results_morpho_am, "Results/phylolm/phylolm_morpho_AM.csv")

head(mat.bm.ea)
dim(mat.bm.ea)
clim.vars <- seq(2,216,2)

mat.bm.ea2 <- mat.bm.ea[clim.vars,]
dim(mat.bm.ea2)
mat.ou.ea2 <- mat.ou.ea[clim.vars,]
dim(mat.ou.ea2)

# write.csv(mat.bm.ea2, "Results/phylolm/coeff_BM_morpho_EA.csv")
# write.csv(mat.ou.ea2, "Results/phylolm/coeff_OU_morpho_EA.csv")


mat.bm.am2 <- mat.bm.am[clim.vars,]
dim(mat.bm.am2)
mat.ou.am2 <- mat.ou.am[clim.vars,]
dim(mat.ou.am2)
# write.csv(mat.bm.am2, "Results/phylolm/coeff_BM_morpho_AM.csv")
# write.csv(mat.ou.am2, "Results/phylolm/coeff_OU_morpho_AM.csv")



# PLOT morpho ~ climate 
##########################
# Plot coefficients
col<- colorRampPalette(c("tomato", "white", "steelblue3"))(100)
heatmap.2(t(mat.am), scale = "none", col = col, 
          trace = "none", density.info = "none", dendrogram= "none", 
          cexCol = 1, Colv=FALSE, Rowv = FALSE, main="American")
pheatmap(t(mat.am), cluster_cols = FALSE,
         cluster_rows = FALSE, border_color = "black", legend = F,
         color = brewer.pal(n = 11, name = "RdBu"))



