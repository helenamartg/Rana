#####################
#     PHYLOLM
# Microhabitat use
#####################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

library(phylolm)
library(phytools)

# 1. Import traits
#===================
dt <- read.csv("Data/dt_59sp.csv", sep=",")
rownames(dt) <- dt$Species

dtEA <- dt[dt$real_reg=="EA", c(2,7,9)]
dtAM <- dt[dt$real_reg=="AM", c(2,7,9)]



# 2. Import phylogeny
#======================
phyloEA <- read.tree("Data/phyloEA_26sp.tre")
phyloAM <- read.tree("Data/phyloAM_33sp.tre")

phyloEA$root.edge <- NULL
phyloAM$root.edge <- NULL


# 3. Check order
#=================
dtEA$Species==phyloEA$tip.label
dtAM$Species==phyloAM$tip.label

dtEA <- dtEA[phyloEA$tip.label,]
dtAM <- dtAM[phyloAM$tip.label,]


# 4. PHYLOLM: SVL ~ microhabitat use
#=====================================
# For Eurasian frogs
model_eurasia_BM <- phylolm(log_svl ~ ECOTYPE, data = dtEA, phy = phyloEA, model = "BM")
model_eurasia_OU <- phylolm(log_svl ~ ECOTYPE, data = dtEA, phy = phyloEA, model = "OUfixedRoot")

# For American frogs
model_america_BM <- phylolm(log_svl ~ ECOTYPE, data = dtAM, phy = phyloAM, model = "BM")
model_america_OU <- phylolm(log_svl ~ ECOTYPE, data = dtAM, phy = phyloAM, model = "OUfixedRoot")

# Compare AICc between models 
AICc_eurasia_BM <- AIC(model_eurasia_BM, corrected=TRUE)
AICc_eurasia_OU <- AIC(model_eurasia_OU, corrected=TRUE)

AICc_america_BM <- AIC(model_america_BM, corrected=TRUE)
AICc_america_OU <- AIC(model_america_OU, corrected=TRUE)

# Choose the best model based on the lower AICc for each radiation
best_model_eurasia <- ifelse(AICc_eurasia_OU < AICc_eurasia_BM, "OU model", "BM model")
best_model_america <- ifelse(AICc_america_OU < AICc_america_BM, "OU model", "BM model")
 
# For both radiations: OU model is the best one

summary(model_eurasia_OU)
summary(model_america_OU)

table(dtEA$ECOTYPE)
table(dtAM$ECOTYPE)


# Calculate deltaAICc
Delta_AIC_BM_Eurasia <- AICc_eurasia_BM - pmin(AICc_eurasia_BM, AICc_eurasia_OU)
Delta_AIC_OU_Eurasia <- AICc_eurasia_OU - pmin(AICc_eurasia_BM, AICc_eurasia_OU)

Delta_AIC_BM_America <- AICc_america_BM - pmin(AICc_america_BM, AICc_america_OU)
Delta_AIC_OU_America <- AICc_america_OU - pmin(AICc_america_BM, AICc_america_OU)


# All results together
results_Eurasia <- cbind(Radiation = "Eurasia", trait = "SVL", AICc_BM = AICc_eurasia_BM, AICc_OU = AICc_eurasia_OU,
                        Alpha = model_eurasia_OU$optpar, DeltaAICc_BM = Delta_AIC_BM_Eurasia, 
                        DeltaAICc_OU = Delta_AIC_OU_Eurasia, Lowest_AICc=NA)

results_America <- cbind(Radiation = "America", trait = "SVL", AICc_BM = AICc_america_BM, AICc_OU = AICc_america_OU,
                         Alpha = model_america_OU$optpar, DeltaAICc_BM = Delta_AIC_BM_America, 
                         DeltaAICc_OU = Delta_AIC_OU_America, Lowest_AICc=NA)

results <- rbind(results_Eurasia, results_America)
results <- as.data.frame(results)

# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results)){
  if (results[i,]$AICc_BM < results[i,]$AICc_OU) {
    results[i,]$Lowest_AICc <- "BM model"
  } else {
    results[i,]$Lowest_AICc <- "OU model"
  }
}


library(tidyr)
library(dplyr)

results[,3] <- as.numeric(results[,3])
results[,4] <- as.numeric(results[,4])
results[,5] <- as.numeric(results[,5])
results[,6] <- as.numeric(results[,6])
results[,7] <- as.numeric(results[,7])
results$Difference_4_units <- abs(results$DeltaAICc_BM - results$DeltaAICc_OU) >= 4


# Select the best model 
results$Best_model <- NA
for (i in 1:nrow(results)){
  if (results[i,9]==TRUE){
    results[i,10] <- results[i,8]
  } else {
    results[i,10] <- "BM model"
  }
}

# coefficients of regressions
bm_Eurasia <- as.data.frame(cbind(summary(model_eurasia_BM)$coefficients, Rsq = model_eurasia_BM$r.squared, 
                          Rsq_adj = model_eurasia_BM$adj.r.squared, 
                          logLik = model_eurasia_BM$logLik, trait="SVL"))

ou_America <- as.data.frame(cbind(summary(model_america_OU)$coefficients, Rsq = model_america_OU$r.squared, 
                                  Rsq_adj = model_america_OU$adj.r.squared, 
                                  logLik = model_america_OU$logLik, trait="SVL"))




# 5. LOOP for all morphometrics
#================================
morpho <- read.csv("Data/morpho_phylores_45sp.csv")
morpho <- morpho[,-1]
morpho <- morpho[morpho$Species!="Rana_tlaloci",]
morpho <- morpho[morpho$Species!="Rana_luteiventris",]
morpho <- morpho[morpho$Species!="Rana_boylii",]
morpho <- morpho[morpho$Species!="Rana_sierrae",]
morpho <- morpho[morpho$Species!="Rana_muscosa",]
morpho <- morpho[morpho$Species!="Rana_aurora",]
morpho <- morpho[morpho$Species!="Rana_cascadae",]
rownames(morpho) <- morpho$Species


# Separate by radiation
########################
library(castor)

# EURASIA
morpho.ea <- morpho[morpho$phylo_reg=="EA",]
morpho.ea <- morpho.ea[,c(1,5,c(12:15))]
results_morpho_ea <- data.frame()
mat.bm.ea <- data.frame()
mat.ou.ea <- data.frame()

for(i in 3:ncol(morpho.ea)) {
  morpho.x <- morpho.ea[complete.cases(morpho.ea[,i]), c(1,2,i)]          # Remove sp from morpho dataset
  phylo.x <- (get_subtree_with_tips(phyloEA, morpho.x$Species))$subtree      # Clip phylogeny
  
  # Set the root edge to 0 (or remove it entirely)
  if(!is.null(phylo.x$root.edge)) {
    phylo.x$root.edge <- NULL  # Remove root edge if it exists
  }
  
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  
  trait <- setNames(morpho.x[,3], morpho.x$Species)
  ecotype <- setNames(morpho.x[,2], morpho.x$Species)
 
    bm_model_eurasia <- phylolm(trait ~ ecotype, phy = phylo.x, model = "BM")
    ou_model_eurasia <- phylolm(trait ~ ecotype, phy = phylo.x, model = "OUfixedRoot")
    
    bm_summary <- summary(bm_model_eurasia)
    ou_summary <- summary(ou_model_eurasia)
    
    bm <- as.data.frame(cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                              Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik))
    bm$trait <- as.character(colnames(morpho.x)[3])
    
    ou <- as.data.frame(cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                              Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik))
    ou$trait <- as.character(colnames(morpho.x)[3])
    
    mat.bm.ea <- rbind(mat.bm.ea, bm)
    mat.ou.ea <- rbind(mat.ou.ea, ou)
    
    AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
    AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)
    
    x <- cbind(Radiation = "Eurasia", 
               trait = colnames(morpho.x)[3], AICc_BM = as.numeric(AICc_bm), AICc_OU = as.numeric(AICc_ou),
               Alpha = ou_summary$optpar, Lowest_AICc=NA)
    
    results_morpho_ea <- rbind(results_morpho_ea, x)
}


# AMERICA
morpho.am <- morpho[morpho$phylo_reg=="AM",]
rownames(morpho.am) <- morpho.am$Species
morpho.am <- morpho.am[,c(1,5,c(12:15))]

results_morpho_am <- data.frame()
mat.bm.am <- data.frame()
mat.ou.am <- data.frame()

for(i in 3:ncol(morpho.am)) {
  morpho.x <- morpho.am[complete.cases(morpho.am[,i]), c(1,2,i)]          # Remove sp from morpho dataset
  phylo.x <- (get_subtree_with_tips(phyloAM, morpho.x$Species))$subtree      # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  
  trait <- setNames(morpho.x[,3], morpho.x$Species)
  ecotype <- setNames(morpho.x[,2], morpho.x$Species)
  
  # Set the root edge to 0 (or remove it entirely)
  if(!is.null(phylo.x$root.edge)) {
    phylo.x$root.edge <- NULL  # Remove root edge if it exists
  }
  
  bm_model_america <- phylolm(trait ~ ecotype, phy = phylo.x, model = "BM")
  ou_model_america <- phylolm(trait ~ ecotype, phy = phylo.x, model = "OUfixedRoot")
  
  bm_summary <- summary(bm_model_america)
  ou_summary <- summary(ou_model_america)
  
  bm <- as.data.frame(cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik))
  bm$trait <- as.character(colnames(morpho.x)[3])
  
  ou <- as.data.frame(cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik))
  ou$trait <- as.character(colnames(morpho.x)[3])
  
  mat.bm.am <- rbind(mat.bm.am, bm)
  mat.ou.am <- rbind(mat.ou.am, ou)
  
  AICc_bm <- AIC(bm_model_america, corrected = TRUE)
  AICc_ou <- AIC(ou_model_america, corrected = TRUE)
  
  x <- cbind(Radiation = "America", 
             trait = colnames(morpho.x)[3], AICc_BM = as.numeric(AICc_bm), AICc_OU = as.numeric(AICc_ou),
             Alpha = ou_summary$optpar, Lowest_AICc=NA)
  
  results_morpho_am <- rbind(results_morpho_am, x)
}


# combine
results_morpho <- rbind(results_morpho_ea, results_morpho_am)
results_morpho <- as.data.frame(results_morpho)
results_morpho[,3] <- as.numeric(results_morpho[,3])
results_morpho[,4] <- as.numeric(results_morpho[,4])
results_morpho[,5] <- as.numeric(results_morpho[,5])

# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results_morpho)){
  if (results_morpho[i,]$AICc_BM < results_morpho[i,]$AICc_OU) {
    results_morpho[i,]$Lowest_AICc <- "BM model"
  } else {
    results_morpho[i,]$Lowest_AICc <- "OU model"
  }
}

# Calculate delta AICc
results_morpho <- results_morpho %>%
  mutate(
    Delta_AIC_BM = AICc_BM - pmin(AICc_BM, AICc_OU),
    Delta_AIC_OU = AICc_OU - pmin(AICc_BM, AICc_OU)
  )

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is  4 units
results_morpho <- results_morpho %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)


# Select the best model 
results_morpho$Best_model <- NA
for (i in 1:nrow(results_morpho)){
  if (results_morpho[i,9]==TRUE){
    results_morpho[i,10] <- results_morpho[i,6]
  } else {
    results_morpho[i,10] <- "BM model"
  }
}

results_morpho
results <- results[,c(1:5,8,6,7,9,10)]

colnames(results) <- colnames(results_morpho)
all.results <- rbind(results, results_morpho)


bm_Eurasia
ou_America


# SAVE
# write.csv(all.results, "Results/phylolm/Ecotype/Ecotype_phylolm_models.csv")
# write.csv(mat.bm.ea, "Results/phylolm/Ecotype/Ecotype_coef_phylolm_morpho_BM_Eurasia.csv")
# write.csv(mat.ou.ea, "Results/phylolm/Ecotype/Ecotype_coef_phylolm_morpho_OU_Eurasia.csv")
# write.csv(mat.ou.am, "Results/phylolm/Ecotype/Ecotype_coef_phylolm_morpho_OU_America.csv")
# write.csv(mat.bm.am, "Results/phylolm/Ecotype/Ecotype_coef_phylolm_morpho_BM_America.csv")
# 
# write.csv(bm_Eurasia, "Results/phylolm/Ecotype/Ecotype_coef_phylolm_svl_BM_Eurasia.csv")
# write.csv(ou_America, "Results/phylolm/Ecotype/Ecotype_coef_phylolm_svl_OU_America.csv")
