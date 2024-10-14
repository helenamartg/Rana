#################################
# Phylogenetic linear models
# with BM and OU: SVL
#################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(phylolm)
library(phytools)


# 1. Import traits
#===================
dt <- read.csv("DATA/dt_68sp.csv", sep=",")
rownames(dt) <- dt$sp

svlEA <- dt[dt$phylo_reg=="EA", c(2,10)]
svlAM <- dt[dt$phylo_reg=="AM", c(2,10)]


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
svlEA <- svlEA[phyloEA$tip.label,]
svlAM <- svlAM[phyloAM$tip.label,]

svlEA$sp==phyloEA$tip.label
svlAM$sp==phyloAM$tip.label

rownames(climEA)==phyloEA$tip.label
rownames(climAM)==phyloAM$tip.label


# 5. PHYLOLM(): Fits a phylogenetic linear regression model. 
#===============

# 5.1 Eurasia
climEA <- as.data.frame(climEA)
svlEA <- setNames(svlEA$log_svl, svlEA$sp)

dt_eurasia <- cbind(svlEA, climEA)
results <- data.frame()
mat.bm <- data.frame()
mat.ou <- data.frame()

for(i in 1:ncol(climEA)) {
    bm_model_eurasia <- phylolm(svlEA ~ climEA[,i], phy = phyloEA, model = "BM")
    ou_model_eurasia <- phylolm(svlEA ~ climEA[,i], phy = phyloEA, model = "OUfixedRoot")
    
    bm_summary <- summary(bm_model_eurasia)
    ou_summary <- summary(ou_model_eurasia)
    
    bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                    Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)
    rownames(bm)[2] <- names(climEA[i])
    
    
    ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                    Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)
    rownames(ou)[2] <- names(climEA[i])
    
    mat.bm <- rbind(mat.bm, bm)
    mat.ou <- rbind(mat.ou, ou)
    
    AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
    AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)
    
    results <- rbind(results, data.frame(names(climEA[i]), AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                         Radiation = "Eurasian", Alpha = ou_summary$optpar, Best_model=NA))
}


# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results)){
  if (results[i,]$AICc_BM < results[i,]$AICc_OU) {
    results[i,]$Best_model <- "BM model"
  } else {
    results[i,]$Best_model <- "OU model"
  }
}


# Calculate delta AICc
library(tidyr)
library(dplyr)
aicc_data <- results %>%
  mutate(
    Delta_AIC_BM = AICc_BM - pmin(AICc_BM, AICc_OU),
    Delta_AIC_OU = AICc_OU - pmin(AICc_BM, AICc_OU)
  )

# View the updated data frame with ΔAIC values
print(aicc_data)

# If deltaAICc > 2, you can confidently say that the model with the lower AICc is the better fit

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)

# Despite OU models have lower AICc, there are no significant differences from BM models 
# so we choose the simplest ones which are BM in order to avoid overfit the data with OU in this case

########
# SAVE 
#########
clim.vars <- seq(2,54,2)
mat.bm <- mat.bm[clim.vars,]
# write.csv(aicc_data, "Results/phylolm/phylolm_SVL_Eurasia.csv")
# write.csv(mat.bm, "Results/phylolm/output_BM_Eurasia.csv")
# write.csv(mat.ou, "Results/phylolm/output_OU_Eurasia.csv")


# 5.2 America
climAM <- as.data.frame(climAM)
svlAM <- setNames(svlAM$log_svl, svlAM$sp)

dt_america <- cbind(svlAM, climAM)
results2 <- data.frame()
mat.bm2 <- data.frame()
mat.ou2 <- data.frame()

for(i in 1:ncol(climAM)) {
  bm_model_america <- phylolm(svlAM ~ climAM[,i], phy = phyloAM, model = "BM")
  ou_model_america <- phylolm(svlAM ~ climAM[,i], phy = phyloAM, model = "OUfixedRoot")
  
  bm_summary <- summary(bm_model_america)
  ou_summary <- summary(ou_model_america)
  
  bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
              Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)
  rownames(bm)[2] <- names(climAM[i])
  
  ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
              Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)
  rownames(ou)[2] <- names(climAM[i])
  
  mat.bm2 <- rbind(mat.bm2, bm)
  mat.ou2 <- rbind(mat.ou2, ou)
  
  AICc_bm <- AIC(bm_model_america, corrected = TRUE)
  AICc_ou <- AIC(ou_model_america, corrected = TRUE)
  
  results2 <- rbind(results2, data.frame(names(climAM[i]), AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                         Radiation = "American", Alpha = ou_summary$optpar, Best_model=NA))
  
}

# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results2)){
  if (results2[i,]$AICc_BM < results2[i,]$AICc_OU) {
    results2[i,]$Best_model <- "BM model"
  } else {
    results2[i,]$Best_model <- "OU model"
  }
}

# Calculate delta AICc
aicc_data2 <- results2 %>%
  mutate(
    Delta_AIC_BM = AICc_BM - pmin(AICc_BM, AICc_OU),
    Delta_AIC_OU = AICc_OU - pmin(AICc_BM, AICc_OU)
  )

# View the updated data frame with ΔAIC values
print(aicc_data2)

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data2 <- aicc_data2 %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data2 <- aicc_data2 %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)


# View the updated data frame with the flag
print(aicc_data2)

# Because deltaAICc > 2 in all cases, the model with the lower AICc (OU) is the better fit.

########
# SAVE 
#########
clim.vars
mat.ou2 <- mat.ou2[clim.vars,]
# write.csv(aicc_data2, "Results/phylolm/aphylo_SVL_America.csv")
# write.csv(mat.bm2, "Results/phylolm/coeff_BM_America.csv")
# write.csv(mat.ou2, "Results/phylolm/coeff_OU_America.csv")


# PLOT RESULTS
#===============
# Eurasia --> BM
mat.bm
eurasia <- mat.bm[,c(1,4)]
eurasia <- cbind(rownames(eurasia), eurasia)
names(eurasia) <- c("clim.var", "svl", "p-value")

coef.mat <- matrix(nrow=nrow(eurasia), ncol = 2)

for( i in 1:nrow(eurasia)) {
  coef.mat[i,1] <- eurasia[i,1]
  if (eurasia$`p-value`[i] > 0.05) {
    coef.mat[i,2] <- 0
  } else {
    coef.mat[i,2] <- round((eurasia$svl[i]),3)
  }
  
}

coef.mat <- as.data.frame(coef.mat)
colnames(coef.mat) <- c("clim.var", "svl")
rownames(coef.mat) <- coef.mat$clim.var
coef.mat$svl <- as.numeric(coef.mat$svl)

library(ggplot2)
ggplot(coef.mat, aes(x=clim.var, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + labs(x= "", y="coefficient") +
  theme_classic() + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("Eurasia") +
  geom_hline(yintercept = 0, color="gray50")+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=15))




# America --> OU
america <- mat.ou2[,c(1,4)]
america <- cbind(rownames(america), america)
names(america) <- c("clim.var", "svl", "p-value")

coef.mat <- matrix(nrow=nrow(america), ncol = 2)

for( i in 1:nrow(america)) {
  coef.mat[i,1] <- america[i,1]
  if (america$`p-value`[i] > 0.05) {
    coef.mat[i,2] <- 0
  } else {
    coef.mat[i,2] <- round((america$svl[i]),3)
  }
  
}

coef.mat <- as.data.frame(coef.mat)
colnames(coef.mat) <- c("clim.var", "svl")
rownames(coef.mat) <- coef.mat$clim.var
coef.mat$svl <- as.numeric(coef.mat$svl)

ggplot(coef.mat, aes(x=clim.var, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + labs(x= "", y="coefficient") +
  theme_classic() + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("America") +
  geom_hline(yintercept = 0, color="gray50")+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=15))


# Correlation between climatic vars
library(caret)
cor_matrix <- cor(climEA)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8)
red_climEA <- climEA[, -high_corr]

library(caret)
cor_matrix <- cor(climAM)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8)
red_climAM <- climAM[, -high_corr]

