#####################################
# Phylogenetic linear models with
# PHYLOLM: svl ~ climatic variables
####################################

# For these analyses I used max SVL data. 

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

library(phylolm)
library(phytools)


# 1. Import traits
#===================
dt <- read.csv("Data/dt_59sp.csv", sep=",")
rownames(dt) <- dt$Species

svlEA <- dt[dt$real_reg=="EA", c(2,9)]
svlAM <- dt[dt$real_reg=="AM", c(2,9)]


# 2. Import phylogeny
#======================
phyloEA <- read.tree("Data/phyloEA_26sp.tre")
phyloAM <- read.tree("Data/phyloAM_33sp.tre")


# 3. Import climatic variables
#===============================
climEA <- read.csv("Data/ClimEA_26sp.csv")
climAM <- read.csv("Data/ClimAM_33sp.csv")

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

svlEA$Species==phyloEA$tip.label
svlAM$Species==phyloAM$tip.label

rownames(climEA)==phyloEA$tip.label
rownames(climAM)==phyloAM$tip.label

climEA <- climEA[match(phyloEA$tip.label, rownames(climEA)),]
climAM <- climAM[phyloAM$tip.label,]



# 5. PHYLOLM(): Fits a phylogenetic linear regression model 
#=============================================================
#######################
# 5.1 Bergmann's rule:   SVL ~ TAVGavg + SRADmin + PETmax
######################

# 1). Eurasia
#=============
climEA1 <- as.data.frame(climEA)
svlEA <- setNames(svlEA$log_svl, svlEA$Species)
dt_eurasia <- cbind(svlEA, climEA1)

results <- data.frame()

phyloEA$root.edge <- NULL
bm_model_eurasia <- phylolm(svlEA ~ TMINavg + SRADavg + PETmax, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "BM")
ou_model_eurasia <- phylolm(svlEA ~ TMINavg + SRADavg + PETmax, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "OUfixedRoot")

bm_summary <- summary(bm_model_eurasia)
ou_summary <- summary(ou_model_eurasia)

  
bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
              Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
              Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)
  
results <- rbind(results, data.frame(model= "bergmann's rule", vars="TMINavg + SRADavg + PETmax",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                       Radiation = "Eurasian", Alpha = ou_summary$optpar, Best_model=NA))

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

# There are no significant differences from BM models 
# so we choose the simplest ones, which are BM in order to avoid overfit the data with OU in this case


########
# SAVE 
#########
BM_bergmann_Eurasia <- bm_model_eurasia
OU_bergmann_Eurasia <- ou_model_eurasia
write.csv(bm, "Results/phylolm/output_BM_bergmann_Eurasia.csv")
write.csv(ou, "Results/phylolm/output_OU_bergmann_Eurasia.csv")




# 2). America
#=============
climAM1 <- as.data.frame(climAM)
svlAM <- setNames(svlAM$log_svl, svlAM$Species)
dt_america <- cbind(svlAM, climAM1)


phyloAM$root.edge <- NULL
bm_model_america <- phylolm(svlAM ~ TMINavg + SRADavg + PETmax, data=dt_america,
                            phy = phyloAM, model = "BM")
ou_model_america <- phylolm(svlAM ~ TMINavg + SRADavg + PETmax, data=dt_america,
                            phy = phyloAM, model = "OUfixedRoot")

bm_summary <- summary(bm_model_america)
ou_summary <- summary(ou_model_america)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_america, corrected = TRUE)
AICc_ou <- AIC(ou_model_america, corrected = TRUE)

results <- rbind(results, data.frame(model= "bergmann's rule", vars="TMINavg + SRADavg + PETmax", 
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "American", Alpha = ou_summary$optpar, Best_model=NA))

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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)


########
# SAVE 
########
BM_bergmann_America <- bm_model_america
OU_bergmann_America <- ou_model_america
write.csv(bm, "Results/phylolm/output_BM_bergmann_America.csv")
write.csv(ou, "Results/phylolm/output_OU_bergmann_America.csv")



####################################
# 5.2 Water conservation hypothesis:   SVL ~ ARIDmin + ARIDmax + PRECmax + PRECmin
####################################         

# 1). Eurasia
#=============
climEA2 <- as.data.frame(climEA)
svlEA
dt_eurasia <- cbind(svlEA, climEA2)
dt_eurasia <- as.data.frame(dt_eurasia)

bm_model_eurasia <- phylolm(svlEA ~ ARIDmax + ARIDmin +
                              PRECmax + PRECmin, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "BM")

ou_model_eurasia <- phylolm(svlEA ~ ARIDmax + ARIDmin +
                              PRECmax + PRECmin, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "OUfixedRoot")

bm_summary <- summary(bm_model_eurasia)
ou_summary <- summary(ou_model_eurasia)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)

results <- rbind(results, data.frame(model= "water conservation hypothesis", 
                                     vars="ARIDmax + ARIDmin + PRECmax + PRECmin",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "Eurasian", Alpha = ou_summary$optpar, Best_model=NA))

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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)


# SAVE
BM_waterconservation_Eurasia <- bm_model_eurasia
OU_waterconservation_Eurasia <- ou_model_eurasia
write.csv(bm, "Results/phylolm/output_BM_waterconservation_Eurasia.csv")
write.csv(ou, "Results/phylolm/output_OU_waterconservation_Eurasia.csv")



# 2). America
#=============
climAM2 <- climAM[,c(1:3, 9:12, 20:21)]
svlAM
dt_america <- cbind(svlAM, climAM2)
dt_america <- as.data.frame(dt_america)

bm_model_america <- phylolm(svlAM ~ ARIDmax + ARIDmin +
                              PRECmax + PRECmin, 
                            data=dt_america,
                            phy = phyloAM, model = "BM")

ou_model_america <- phylolm(svlAM ~ ARIDmax + ARIDmin +
                              PRECmax + PRECmin, 
                            data=dt_america,
                            phy = phyloAM, model = "OUfixedRoot")

bm_summary <- summary(bm_model_america)
ou_summary <- summary(ou_model_america)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_america, corrected = TRUE)
AICc_ou <- AIC(ou_model_america, corrected = TRUE)

results <- rbind(results, data.frame(model= "water conservation hypothesis", 
                                     vars="ARIDmax + ARIDmin + PRECmax + PRECmin",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "America", Alpha = ou_summary$optpar, Best_model=NA))

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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)

# SAVE
BM_waterconservation_America <- bm_model_america
OU_waterconservation_America <- ou_model_america
write.csv(bm, "Results/phylolm/output_BM_waterconservation_America.csv")
write.csv(ou, "Results/phylolm/output_OU_waterconservation_America.csv")



#######################################
# 5.2 Primary productivity hypothesis:   SVL ~ NDVImax + NDVImin
####################################### 

# 1). Eurasia
#==============
climEA3 <- climEA[,4:6]
svlEA
dt_eurasia <- cbind(svlEA, climEA3)
dt_eurasia <- as.data.frame(dt_eurasia)

bm_model_eurasia <- phylolm(svlEA ~ NDVImax + NDVImin, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "BM")

ou_model_eurasia <- phylolm(svlEA ~ NDVImax + NDVImin, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "OUfixedRoot")

bm_summary <- summary(bm_model_eurasia)
ou_summary <- summary(ou_model_eurasia)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)

results <- rbind(results, data.frame(model= "primary productivity hypothesis", 
                                     vars="NDVImax + NDVImin",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "Eurasian", Alpha = ou_summary$optpar, Best_model=NA))

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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)


# SAVE
BM_productivity_Eurasia <- bm_model_eurasia
OU_productivity_Eurasia <- ou_model_eurasia
write.csv(bm, "Results/phylolm/output_BM_productivity_Eurasia.csv")
write.csv(ou, "Results/phylolm/output_OU_productivity_Eurasia.csv")


# 2). America
#==============
climAM3 <- climAM[,4:6]
svlAM
dt_america <- cbind(svlAM, climAM3)
dt_america <- as.data.frame(dt_america)

bm_model_america <- phylolm(svlAM ~ NDVImax + NDVImin, 
                            data=dt_america,
                            phy = phyloAM, model = "BM")

ou_model_america <- phylolm(svlAM ~ NDVImax + NDVImin, 
                            data=dt_america,
                            phy = phyloAM, model = "OUfixedRoot")

bm_summary <- summary(bm_model_america)
ou_summary <- summary(ou_model_america)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_america, corrected = TRUE)
AICc_ou <- AIC(ou_model_america, corrected = TRUE)

results <- rbind(results, data.frame(model= "primary productivity hypothesis", 
                                     vars="NDVImax + NDVImin",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "American", Alpha = ou_summary$optpar, Best_model=NA))

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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)


# SAVE
BM_productivity_America <- bm_model_america
OU_productivity_America <- ou_model_america
write.csv(bm, "Results/phylolm/output_BM_productivity_American.csv")
write.csv(ou, "Results/phylolm/output_OU_productivity_American.csv")

aicc_data <- rbind(aicc_data, aicc_data)


###########################
# 5.4 Habitat availability:   SVL ~ ELEVmax + ELEVmin
###########################

# 1). Eurasia
#==============
climEA4 <- climEA[,23:24]
svlEA
dt_eurasia <- cbind(svlEA, climEA4)
dt_eurasia <- as.data.frame(dt_eurasia)

bm_model_eurasia <- phylolm(svlEA ~ ELEVmax + ELEVmin, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "BM")

ou_model_eurasia <- phylolm(svlEA ~ ELEVmax + ELEVmin, 
                            data=dt_eurasia,
                            phy = phyloEA, model = "OUfixedRoot")

bm_summary <- summary(bm_model_eurasia)
ou_summary <- summary(ou_model_eurasia)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)

results <- rbind(results, data.frame(model= "habitat availability",
                                     vars="ELEVmax + ELEVmin",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "Eurasian", Alpha = ou_summary$optpar, Best_model=NA))


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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)


#SAVE
BM_habitat_Eurasia <- bm_model_eurasia
OU_habitat_Eurasia <- ou_model_eurasia
write.csv(bm, "Results/phylolm/output_BM_habitat_Eurasia.csv")
write.csv(ou, "Results/phylolm/output_OU_habitat_Eurasia.csv")



# 2). America
#==============
climAM4 <- climAM[,23:24]
svlAM
dt_america <- cbind(svlAM, climAM4)
dt_america <- as.data.frame(dt_america)

bm_model_america <- phylolm(svlAM ~ ELEVmax + ELEVmin, 
                            data=dt_america,
                            phy = phyloAM, model = "BM")

ou_model_america <- phylolm(svlAM ~ ELEVmax + ELEVmin, 
                            data=dt_america,
                            phy = phyloAM, model = "OUfixedRoot")

bm_summary <- summary(bm_model_america)
ou_summary <- summary(ou_model_america)


bm <- cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik)

ou <- cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik)

AICc_bm <- AIC(bm_model_america, corrected = TRUE)
AICc_ou <- AIC(ou_model_america, corrected = TRUE)

results <- rbind(results, data.frame(model= "habitat availability",
                                     vars="ELEVmax + ELEVmin",
                                     AICc_BM = AICc_bm, AICc_OU = AICc_ou,
                                     Radiation = "America", Alpha = ou_summary$optpar, Best_model=NA))


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

# Add a new column to check if the difference between Delta_AIC_BM and Delta_AIC_OU is >= 2 units
aicc_data <- aicc_data %>%
  mutate(Difference_2_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 2)

aicc_data <- aicc_data %>%
  mutate(Difference_4_units = abs(Delta_AIC_BM - Delta_AIC_OU) >= 4)

# View the updated data frame with the flag
print(aicc_data)


BM_habitat_America <- bm_model_america
OU_habitat_America <- ou_model_america

# SAVE
write.csv(bm, "Results/phylolm/output_BM_habitat_American.csv")
write.csv(ou, "Results/phylolm/output_OU_habitat_American.csv")

write.csv(aicc_data, "Results/phylolm/output_all_models.csv")




#################################################################################
# PLOT 

library(ggplot2)
library(dplyr)


# Extraer coeficientes y p-valores de los modelos
extraer_coef <- function(modelo, nombre_modelo) {
  data.frame(
    Variable = names(modelo$coefficients)[-1],  # Excluir el intercepto
    Estimate = modelo$coefficients[-1],         # Coeficientes sin el intercepto
    p_value = summary(modelo)$coefficients[-1,4],   # Valores p
    Modelo = nombre_modelo
  )
}

# 2) Modelos EURASIA
#####################
coef_data <- bind_rows(
  extraer_coef(BM_bergmann_Eurasia, "Bergmann's rule"),
  extraer_coef(BM_waterconservation_Eurasia, "Water conservation hypothesis"),
  extraer_coef(OU_productivity_Eurasia, "Primary productivity hypothesis"),
  extraer_coef(BM_habitat_Eurasia, "Habitat availability hypothesis")
) %>%
  mutate(Significativo = ifelse(p_value < 0.05, "Sí", "No"))

coef_data$Estimate <- ifelse(coef_data$p_value > 0.05, 0, round(coef_data$Estimate, 3))

# Gráfico de coeficientes
coef_data$Variable <- factor(coef_data$Variable, levels = c(
  "PETmax", "SRADavg", "TMINavg", # Bergmann's rule
  "ARIDmax", "ARIDmin", "PRECmax", "PRECmin", # Water conservation hypothesis
  "NDVImax", "NDVImin", # Primary productivity hypothesis
  "ELEVmax", "ELEVmin"  # Habitat availability hypothesis
))
p <- ggplot(coef_data, aes(x=Variable, y=Estimate, fill=Estimate > 0)) +
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



# 2) Modelos AMERICA
#####################
coef_data <- bind_rows(
  extraer_coef(OU_bergmann_America, "Bergmann's rule"),
  extraer_coef(OU_waterconservation_America, "Water conservation hypothesis"),
  extraer_coef(OU_productivity_America, "Primary productivity hypothesis"),
  extraer_coef(OU_habitat_America, "Habitat availability hypothesis")
) %>%
  mutate(Significativo = ifelse(p_value < 0.05, "Sí", "No"))

coef_data$Estimate <- ifelse(coef_data$p_value > 0.05, 0, round(coef_data$Estimate, 3))

# Gráfico de coeficientes
coef_data$Variable <- factor(coef_data$Variable, levels = c(
  "PETmax", "SRADavg", "TMINavg", # Bergmann's rule
  "ARIDmax", "ARIDmin", "PRECmax", "PRECmin", # Water conservation hypothesis
  "NDVImax", "NDVImin", # Primary productivity hypothesis
  "ELEVmax", "ELEVmin"  # Habitat availability hypothesis
))

q <- ggplot(coef_data, aes(x=Variable, y=Estimate, fill=Estimate > 0)) +
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


# Arrange in two columns
library(gridExtra)
grid.arrange(p, q, ncol=1)


