##############################
# Multivariate PGLS
#############################

rm(list=ls())
library(phytools)
library(ape)
library(geomorph)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")


# 1. Import phylogeny
#=======================
phyloEA <- read.tree("DATA/phyloEA_31sp.tre")
phyloAM <- read.tree("DATA/phyloAM_37sp.tre")


# 2. Import trait database
#======================================
dt68 <- read.csv("DATA/dt_68sp.csv")
dt68 <- dt68[,-1]
rownames(dt68) <- dt68$sp

# Split dataset
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

# Same order as phylogeny
dtEA <- dtEA[phyloEA$tip.label,]
dtAM <- dtAM[phyloAM$tip.label,]


# 3. Import climatic data
#=========================
climEA <- read.csv("DATA/clim_EA_31sp.csv")
climAM <- read.csv("DATA/clim_AM_37sp.csv")

rownames(climEA) <- climEA$X
climEA <- climEA[,-c(1,2)]

rownames(climAM) <- climAM$X
climAM <- climAM[,-c(1,2)]

# Remove SD and scale
climEA.s <- scale(climEA[,-grep("sd", colnames(climEA))])
climAM.s <- scale(climAM[,-grep("sd", colnames(climAM))])

# MAtch phylogeny
climEA.s <- climEA.s[phyloEA$tip.label,]
climAM.s <- climAM.s[phyloAM$tip.label,]


# 4. Select significant variables only 
#======================================
coef.ea <- read.csv("Results/PGLS/Eurasia/coef_matrix_EA.txt")
ncol(coef.ea)==ncol(climEA.s)

sign.clim.ea <- NULL
for (i in 1:ncol(coef.ea)){
  if (coef.ea[1,i]!=0){
    x <- as.matrix(climEA.s[,i])
    colnames(x) <- colnames(coef.ea[i])
    sign.clim.ea <- cbind(sign.clim.ea,x)
  }
}

coef.am <- read.csv("Results/PGLS/America/coef_matrix_AM.txt")
ncol(coef.am)==ncol(climAM.s)

sign.clim.am <- NULL
for (i in 1:ncol(coef.am)){
  if (coef.am[1,i]!=0){
    x <- as.matrix(climAM.s[,i])
    colnames(x) <- colnames(coef.am[i])
    sign.clim.am <- cbind(sign.clim.am,x)
  }
}

# Remove repeated vars
# SRADmax, PETmax, NDVImax, ELEVmax, ARIDmax
sign.clim.am <- sign.clim.am[,-c(1,4,6,10,14)]


# 5. PGLS: lm.rrpp {SVL ~ significant clim vars + ecotype + biome}
#===================================================================
# Build the Phylogenetic var-covariance matrix (correlation)
cov.ea <- vcv(phyloEA, model="Brownian") 
cov.am <- vcv(phyloAM, model="Brownian") 


#####################################
# 5.1 INDIVIDUAL CLIM EURASIA
#####################################
dt.clim <- rrpp.data.frame(sign.clim.ea)
dt.all <- rrpp.data.frame(clim = dt.clim, svl = log(dtEA$SVL), biome = dtEA$biome, 
                          ecotype = dtEA$ecotype, phy.cvc = cov.ea, clutch=log(dtEA$clutch_size))
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(sign.clim.ea)]))

# model : SVL ~ clim vars + biome + eoctype + clutch

for (i in 1:length(dt.all[1:ncol(sign.clim.ea)])){
  
  dt.temp <- rrpp.data.frame(x1 = dt.all[[1]][,i], 
                             biome = dt.all$biome,
                             ecotype = dt.all$ecotype,
                             clutch = dt.all$clutch,
                             y = dt.all$svl, phy.cvc = cov.ea)
  
  fitOLS <- lm.rrpp(y ~ x1 + biome + ecotype + clutch, data = dt.temp, print.progress = FALSE, iter=999)
  fitGLS <- lm.rrpp(y ~ x1 + biome + ecotype + clutch, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
  ols <- cbind(anova(fitOLS)$table, colnames(sign.clim.ea)[i],"svl")
  gls <- cbind(anova(fitGLS)$table, colnames(sign.clim.ea)[i],"svl")
  
  anovas_OLS <- rbind(anovas_OLS, ols)
  anovas_GLS <- rbind(anovas_GLS, gls)
  
  coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(sign.clim.ea)[i], 
                                "biome", "ecotype","svl"))
  
  if (anova(fitGLS)$table[1,7] > 0.05) {
    coef.matrix[1,i] <- 0
  } else {
    coef.matrix[1,i] <- (fitGLS$LM$gls.coefficients[2])
  }
  
}
colnames(coef.matrix) <- c(colnames(sign.clim.ea))
rownames(coef.matrix) <- "svl"
coef.matrix <- round(coef.matrix, 3)
anovas_GLS
anovas_OLS
colnames(coef_GLS) <- c("intercept", "x", "clim_var", "biome", "ecotype", "svl")
# No significant relationships
# write.csv(anovas_GLS, "anovas_GLS_mult_EA.csv")


###############################
# 5.1 INDIVIDUAL CLIM AMERICA
################################
dt.clim <- rrpp.data.frame(sign.clim.am)
dt.all <- rrpp.data.frame(clim = dt.clim, svl = log(dtAM$SVL), biome = dtAM$biome, 
                          ecotype = dtAM$ecotype, clutch= log(dtAM$clutch_size), phy.cvc = cov.am)
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(sign.clim.am)]))

# model : SVL ~ clim vars + biome + eoctype + clutch

for (i in 1:length(dt.all[1:ncol(sign.clim.am)])){
  
  dt.temp <- rrpp.data.frame(x1 = dt.all[[1]][,i], 
                             biome = dt.all$biome,
                             ecotype = dt.all$ecotype,
                             clutch = dt.all$clutch,
                             y = dt.all$svl, phy.cvc = cov.am)
  
  fitOLS <- lm.rrpp(y ~ x1 + biome + ecotype + clutch, data = dt.temp, print.progress = FALSE, iter=999)
  fitGLS <- lm.rrpp(y ~ x1 + biome + ecotype + clutch, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
  ols <- cbind(anova(fitOLS)$table, colnames(sign.clim.am)[i],"svl")
  gls <- cbind(anova(fitGLS)$table, colnames(sign.clim.am)[i],"svl")
  
  anovas_OLS <- rbind(anovas_OLS, ols)
  anovas_GLS <- rbind(anovas_GLS, gls)
  
  coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(sign.clim.am)[i], 
                                "biome", "ecotype","svl"))
  
  if (anova(fitGLS)$table[1,7] > 0.05) {
    coef.matrix[1,i] <- 0
  } else {
    coef.matrix[1,i] <- (fitGLS$LM$gls.coefficients[2])
  }
  
}
colnames(coef.matrix) <- c(colnames(sign.clim.am), "biome", "ecotype")
rownames(coef.matrix) <- "svl"
coef.matrix <- round(coef.matrix, 3)
anovas_GLS
anovas_OLS
colnames(coef_GLS) <- c("intercept", "x", "clim_var", "svl")
# Some significant relationships between some climatic vars + clutch ~ svl 
# write.csv(anovas_GLS, "anovas_GLS_mult_AM.csv")

#=================================
# 6. Models BY HAND without loop
#================================
# EURASIA
#========
# SVL ~ SRADmin + SRADavg + PETmin + PETavg + NDVImin + ARIDmax + biome + ecotype
sign.clim.ea <- as.data.frame(sign.clim.ea)
dt.all <- rrpp.data.frame(SRADmin = sign.clim.ea$SRADmin,
                          SRADavg = sign.clim.ea$SRADavg,
                          PETmin = sign.clim.ea$PETmin,
                          PETavg = sign.clim.ea$PETavg,
                          NDVImin = sign.clim.ea$NDVImin,
                          ARIDmax = sign.clim.ea$ARIDmax,
                          svl = dtEA$SVL, 
                          biome = dtEA$biome, 
                          ecotype = dtEA$ecotype, phy.cvc = cov.ea)
str(dt.all)

# FULL models
fitOLS <- lm.rrpp(svl ~ SRADmin + SRADavg + PETmin + PETavg + NDVImin + ARIDmax 
                  + biome + ecotype, data = dt.all, print.progress = FALSE, iter=999)

fitGLS <- lm.rrpp(svl ~ SRADmin + SRADavg + PETmin + PETavg + NDVImin + ARIDmax 
                  + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

anova(fitOLS)$table
anova(fitGLS)$table


# Separated models (GLS)
fit1 <- lm.rrpp(svl ~ SRADmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit2 <- lm.rrpp(svl ~ SRADavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit3 <- lm.rrpp(svl ~ PETmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit4 <- lm.rrpp(svl ~ PETavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit5 <- lm.rrpp(svl ~ NDVImin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit6 <- lm.rrpp(svl ~ ARIDmax + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit7 <- lm.rrpp(svl ~ SRADmin + SRADavg + PETmin + PETavg + NDVImin + ARIDmax 
                          + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit8 <- lm.rrpp(svl ~ ecotype + biome, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit9 <- lm.rrpp(svl ~ biome, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)


# Compare among models: AIC
y <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9,  type = "logLik", tol = 0.01)
y <- setNames(y$table$AIC, y$names)
sort(y) # Best model: SVL ~ SRADmin + biome + ecotype


# AMERICA
#==========
# SVL ~ WINDmin + WINDmax  + WINDavg + SRADmin + PRECmin + PRECavg + PETmin
#       + NDVImin+ NDVIavg + ELEVmin + ELEVavg + ARIDmin

sign.clim.am <- as.data.frame(sign.clim.am)
dt.all <- rrpp.data.frame(WINDmin = sign.clim.am$WINDmin,
                          WINDmax = sign.clim.am$WINDmax,
                          WINDavg = sign.clim.am$WINDavg,
                          SRADmin = sign.clim.am$SRADmin,
                          PRECmin = sign.clim.am$PRECmin,
                          PRECavg = sign.clim.am$PRECavg,
                          PETmin = sign.clim.am$PETmin,
                          NDVImin = sign.clim.am$NDVImin,
                          NDVIavg = sign.clim.am$NDVIavg,
                          ELEVmin = sign.clim.am$ELEVmin,
                          ELEVavg = sign.clim.am$ELEVavg,
                          ARIDmin = sign.clim.am$ARIDmin,
                          svl = dtAM$SVL, 
                          biome = dtAM$biome, 
                          ecotype = dtAM$ecotype, phy.cvc = cov.am)
str(dt.all)

# Full models
fitOLS <- lm.rrpp(svl ~ WINDmin + WINDmax  + WINDavg + SRADmin + PRECmin + PRECavg + PETmin
                  + NDVImin+ NDVIavg + ELEVmin + ELEVavg + ARIDmin 
                  + biome + ecotype, data = dt.all, print.progress = FALSE, iter=999)

fitGLS <- lm.rrpp(svl ~ WINDmin + WINDmax  + WINDavg + SRADmin + PRECmin + PRECavg + PETmin
                  + NDVImin+ NDVIavg + ELEVmin + ELEVavg + ARIDmin 
                  + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

anova(fitOLS)$table
anova(fitGLS)$table

lm.biome <- lm.rrpp(svl ~ biome, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
anova(lm.biome)$table  #Significant


# Separated models (GLS)
fit1 <- lm.rrpp(svl ~ WINDmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit2 <- lm.rrpp(svl ~ WINDmax + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit3 <- lm.rrpp(svl ~ WINDavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit4 <- lm.rrpp(svl ~ SRADmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit5 <- lm.rrpp(svl ~ PRECmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit6 <- lm.rrpp(svl ~ PRECavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit7 <- lm.rrpp(svl ~ PETmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit8 <- lm.rrpp(svl ~ NDVImin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit9 <- lm.rrpp(svl ~ NDVIavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit10 <- lm.rrpp(svl ~ ELEVmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit11 <- lm.rrpp(svl ~ ELEVavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit12 <- lm.rrpp(svl ~ ARIDmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

fit13 <- lm.rrpp(svl ~ WINDmin + WINDmax  + WINDavg + SRADmin + PRECmin + PRECavg + PETmin
                + NDVImin+ NDVIavg + ELEVmin + ELEVavg + ARIDmin 
                + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

fit14 <- lm.rrpp(svl ~ WINDmin + WINDmax  + WINDavg + SRADmin + PRECmin + PRECavg + PETmin
                 + NDVImin+ NDVIavg + ELEVmin + ELEVavg + ARIDmin, 
                 data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

fit15 <- lm.rrpp(svl ~ ecotype + biome, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit16 <- lm.rrpp(svl ~ biome, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)


# Compare models: AIC
x <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9,
                 fit10, fit11, fit12, fit13, fit14, fit15, fit16, type = "logLik", tol = 0.01)

x <- setNames(x$table$AIC, x$names)
sort(x)    # Lowest AIC value for the full model without biome and ecotype



# 7. BOTH RADIATIONS TOGHETER
#=============================
phyloRana <- read.tree("DATA/phyloRana_74sp.tre")
dt74 <- read.csv("DATA/dt_74sp.csv")
biome <- read.csv("DATA/raw_data/biome_all.csv")
biome <- biome[biome$sp_lin!="Rana_psilonota",]

rownames(biome)<- biome$sp_lin
rownames(dt74) <- dt74$sp

biome <- biome[dt74$sp,]
biome$sp_lin==dt74$sp
dt74$biome <- biome$new_biome

dt74 <- dt74[phyloRana$tip.label,]
dt74$sp == phyloRana$tip.label

clim <- read.csv("DATA/clim_74sp.csv")
clim <- clim[,-grep("sd", colnames(clim))]
rownames(clim) <- clim$X
clim <- clim[,-c(1,2)]
clim <- scale(clim)

clim <- clim[phyloRana$tip.label,]


# Prepare dataset for PGLS
cov.rana = vcv(phyloRana, model="Brownian")

# Select only significant climatic variables (from previous PGLS)
clim <- as.data.frame(clim)
dt.all <- rrpp.data.frame(WINDmin = clim$WINDmin,
                          WINDmax = clim$WINDmax,
                          WINDavg = clim$WINDavg,
                          SRADmin = clim$SRADmin,
                          SRADavg = clim$SRADavg,
                          PRECmin = clim$PRECmin,
                          PETmin = clim$PETmin,
                          NDVImin = clim$NDVImin,
                          ELEVmin = clim$ELEVmin,
                          ELEVavg = clim$ELEVavg,
                          ARIDmin = clim$ARIDmin,
                          svl = dt74$SVL, 
                          biome = dt74$biome, 
                          ecotype = dt74$ecotype, phy.cvc = cov.rana)
str(dt.all)

# Full models

fit1 <- lm.rrpp(svl ~ WINDmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit2 <- lm.rrpp(svl ~ WINDmax + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit3 <- lm.rrpp(svl ~ WINDavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit4 <- lm.rrpp(svl ~ SRADmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit5 <- lm.rrpp(svl ~ SRADavg+ biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit6 <- lm.rrpp(svl ~ PRECmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit7 <- lm.rrpp(svl ~ PETmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit8 <- lm.rrpp(svl ~ NDVImin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit9 <- lm.rrpp(svl ~ ELEVmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit10 <- lm.rrpp(svl ~ ELEVavg + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit11 <- lm.rrpp(svl ~ ARIDmin + biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

fit12 <- lm.rrpp(svl ~ biome + ecotype, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)
fit13 <- lm.rrpp(svl ~ biome, data = dt.all, Cov = dt.all$phy.cvc, print.progress = F)

# Compare models (AIC)
z <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9,
                 fit10, fit11, fit12, fit13, type = "logLik", tol = 0.01)
z <- setNames(z$table$AIC, z$names)
sort(z)  # Best model: SVL ~ SRADmin + biome + ecotype
