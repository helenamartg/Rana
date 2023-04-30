##########################################################################
#                       CLIMATIC PGLS
#Does climatic variables have a significant relationship with SVL?
##########################################################################

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
dt74 <- dt74[,-1]

# 3. Import climatic variables
#=============================
clim <- read.csv("DATA/clim_74sp.csv")
clim <- clim[,-c(1)]
rownames(clim) <- clim$X

# 4. Check order
#===============
phyloRana$tip.label == dt74$sp
dt74 <- dt74[phyloRana$tip.label,]

phyloRana$tip.label == clim$X
clim <- clim[phyloRana$tip.label,]


# 5. Scale climatic variables
#============================
# first remove SD
clim2 <- clim[,-grep("sd", colnames(clim))]
clim3 <- scale(clim2[,-1])


# 6. PGLS: lm.rrpp {SVL ~ all climatic variables}
#=================================================
# Build the Phylogenetic var-covariance matrix (correlation)
cov.mat.rana <- vcv(phyloRana, model="Brownian") 

dt.clim <- rrpp.data.frame(clim3)
dt.all <- rrpp.data.frame(clim = dt.clim, svl = dt74$SVL, phy.cvc = cov.mat.rana)
str(dt.all)


anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(clim3)]))

for (i in 1:length(dt.all[1:ncol(clim3)])){
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,i], y = dt.all$svl, phy.cvc = cov.mat.rana)
    fitOLS <- lm.rrpp(y ~ x, data = dt.temp, print.progress = FALSE, iter=999)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
    ols <- cbind(anova(fitOLS)$table, colnames(clim3)[i], "svl")
    gls <- cbind(anova(fitGLS)$table, colnames(clim3)[i], "svl")
    
    anovas_OLS <- rbind(anovas_OLS, ols)
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim3)[i], 
                                  "svl"))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[1,i] <- 0
    } else {
      coef.matrix[1,i] <- (fitGLS$LM$gls.coefficients[2])
    }
    
  }
colnames(coef.matrix) <- colnames(clim3)
rownames(coef.matrix) <- "svl"
coef.matrix <- round(coef.matrix, 3)
anovas_GLS
anovas_OLS
colnames(coef_GLS) <- c("intercept", "x", "clim_var", "svl")

# SAVE
# write.table(anovas_OLS, "anovas_OLS_rana.txt", sep = ",")
# write.table(anovas_GLS, "anovas_GLS_rana.txt", sep = ",")
# write.table(coef_GLS, "coef_GLS_rana.txt", sep = ",")
# write.table(coef.matrix, "coef_matrix_rana.txt", sep =",")

# PLOT COEFFICIENTS
#===================
par(mar=c(6,6,3,6))
library(ggplot2)
coef <- as.data.frame(t(coef.matrix))
coef$vars <- rownames(coef)
ggplot(coef, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "") +
  theme_classic()


# 7. BY radiation
#================
# Import phylogeny
phyloEA <- read.tree("DATA/phyloEA_31sp.tre")
phyloAM <- read.tree("DATA/phyloAM_37sp.tre")

#Import climatic data
climEA <- read.csv("DATA/clim_EA_31sp.csv")
rownames(climEA) <- climEA$X.1
climEA <- climEA[,-c(1,2)]

climAM <- read.csv("DATA/clim_AM_37sp.csv")
rownames(climAM) <- climAM$X.1
climAM <- climAM[,-c(1,2)]

# Same order in phylogeny
climEA <- climEA[phyloEA$tip.label,]
climAM <- climAM[phyloAM$tip.label,]

# Remove SD and scale
climEA <- climEA[,-grep("sd", colnames(climEA))]
climEA.scaled <- scale(climEA)
climAM <- climAM[,-grep("sd", colnames(climAM))]
climAM.scaled <- scale(climAM)

# Import traits
dt68 <- read.csv("DATA/dt_68sp.csv")

dtEA <- dt68[dt68$phylo_reg=="EA",]
dtEA <- dtEA[,-1]
rownames(dtEA) <- dtEA$sp

dtAM <- dt68[dt68$phylo_reg=="AM",]
dtAM <- dtAM[,-1]
rownames(dtAM) <- dtAM$sp

# Same order as phylogeny
dtEA <- dtEA[phyloEA$tip.label,]
dtAM <- dtAM[phyloAM$tip.label,]


#=============
# 7.1 EURASIA
#=============
cov.mat.ea <- vcv(phyloEA, model="Brownian") 

dt.clim.ea <- rrpp.data.frame(climEA.scaled)
dt.all <- rrpp.data.frame(clim = dt.clim.ea, svl = dtEA$SVL, 
                          phy.cvc = cov.mat.ea)
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(climEA.scaled)]))

for (i in 1:length(dt.all[1:ncol(climEA.scaled)])){
  dt.temp <- rrpp.data.frame(x = dt.all[[1]][,i], y = dt.all$svl, phy.cvc = cov.mat.ea)
  fitOLS <- lm.rrpp(y ~ x, data = dt.temp, print.progress = FALSE, iter=999)
  fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
  ols <- cbind(anova(fitOLS)$table, colnames(climEA.scaled)[i], "svl")
  gls <- cbind(anova(fitGLS)$table, colnames(climEA.scaled)[i], "svl")
  
  anovas_OLS <- rbind(anovas_OLS, ols)
  anovas_GLS <- rbind(anovas_GLS, gls)
  
  coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(climEA.scaled)[i], 
                                "svl"))
  
  if (anova(fitGLS)$table[1,7] > 0.05) {
    coef.matrix[1,i] <- 0
  } else {
    coef.matrix[1,i] <- (fitGLS$LM$gls.coefficients[2])
  }
  
}
colnames(coef.matrix) <- colnames(climEA.scaled)
rownames(coef.matrix) <- "svl"
coef.matrix <- round(coef.matrix, 3)
anovas_GLS
anovas_OLS
colnames(coef_GLS) <- c("intercept", "x", "clim_var", "svl")

# SAVE
# write.table(anovas_OLS, "anovas_OLS_EA.txt", sep = ",")
# write.table(anovas_GLS, "anovas_GLS_EA.txt", sep = ",")
# write.table(coef_GLS, "coef_GLS_EA.txt", sep = ",")
# write.table(coef.matrix, "coef_matrix_EA.txt", sep =",")

# Plot coefficients
#===================
par(mar=c(6,6,3,6))
library(ggplot2)
coef <- as.data.frame(t(coef.matrix))
coef$vars <- rownames(coef)
ggplot(coef, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "") +
  theme_classic()


#=============
# 7.2 AMERICA
#=============
cov.mat.am <- vcv(phyloAM, model="Brownian") 

dt.clim.am <- rrpp.data.frame(climAM.scaled)
dt.all <- rrpp.data.frame(clim = dt.clim.am, svl = dtAM$SVL, phy.cvc = cov.mat.am)
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(climAM.scaled)]))

for (i in 1:length(dt.all[1:ncol(climAM.scaled)])){
  dt.temp <- rrpp.data.frame(x = dt.all[[1]][,i], y = dt.all$svl, phy.cvc = cov.mat.am)
  fitOLS <- lm.rrpp(y ~ x, data = dt.temp, print.progress = FALSE, iter=999)
  fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
  ols <- cbind(anova(fitOLS)$table, colnames(climAM.scaled)[i], "svl")
  gls <- cbind(anova(fitGLS)$table, colnames(climAM.scaled)[i], "svl")
  
  anovas_OLS <- rbind(anovas_OLS, ols)
  anovas_GLS <- rbind(anovas_GLS, gls)
  
  coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(climAM.scaled)[i], 
                                "svl"))
  
  if (anova(fitGLS)$table[1,7] > 0.05) {
    coef.matrix[1,i] <- 0
  } else {
    coef.matrix[1,i] <- (fitGLS$LM$gls.coefficients[2])
  }
  
}
colnames(coef.matrix) <- colnames(climAM.scaled)
rownames(coef.matrix) <- "svl"
coef.matrix <- round(coef.matrix, 3)
anovas_GLS
anovas_OLS
colnames(coef_GLS) <- c("intercept", "x", "clim_var", "svl")

# SAVE
# write.table(anovas_OLS, "anovas_OLS_AM.txt", sep = ",")
# write.table(anovas_GLS, "anovas_GLS_AM.txt", sep = ",")
# write.table(coef_GLS, "coef_GLS_AM.txt", sep = ",")
# write.table(coef.matrix, "coef_matrix_AM.txt", sep =",")

# Plot coefficients
#===================
par(mar=c(6,6,3,6))
library(ggplot2)
coef <- as.data.frame(t(coef.matrix))
coef$vars <- rownames(coef)
ggplot(coef, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "") +
  theme_classic()
