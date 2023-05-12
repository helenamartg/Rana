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
dt.all <- rrpp.data.frame(clim = dt.clim, svl = log(dt74$SVL), phy.cvc = cov.mat.rana)
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
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="log(svl)") +
  theme_classic()

# Cambiar los colores del plot
p <- ggplot(coef, aes(x=vars, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="coeffient") +
  theme_classic()

p+ scale_fill_manual(values = c("#4F6367","#EEF5DB"),
                     labels=c('TRUE'='>0','FALSE'='<0'))+
  labs(fill='svl')


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
dt.all <- rrpp.data.frame(clim = dt.clim.ea, svl = log(dtEA$SVL), 
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
# write.table(anovas_OLS, "Results/PGLS/Eurasia/anovas_OLS_EA.txt", sep = ",")
# write.table(anovas_GLS, "Results/PGLS/Eurasia/anovas_GLS_EA.txt", sep = ",")
# write.table(coef_GLS, "Results/PGLS/Eurasia/coef_GLS_EA.txt", sep = ",")
# write.table(coef.matrix, "Results/PGLS/Eurasia/coef_matrix_EA.txt", sep =",")

# Plot coefficients
#===================
par(mar=c(6,6,3,6))
library(ggplot2)
coef <- as.data.frame(t(coef.matrix))
coef$vars <- rownames(coef)
ggplot(coef, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="log(svl)") +
  theme_classic()


#=============
# 7.2 AMERICA
#=============
cov.mat.am <- vcv(phyloAM, model="Brownian") 

dt.clim.am <- rrpp.data.frame(climAM.scaled)
dt.all <- rrpp.data.frame(clim = dt.clim.am, svl = log(dtAM$SVL), phy.cvc = cov.mat.am)
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
# write.table(anovas_OLS, "Results/PGLS/America/anovas_OLS_AM.txt", sep = ",")
# write.table(anovas_GLS, "Results/PGLS/America/anovas_GLS_AM.txt", sep = ",")
# write.table(coef_GLS, "Results/PGLS/America/coef_GLS_AM.txt", sep = ",")
# write.table(coef.matrix, "Results/PGLS/America/coef_matrix_AM.txt", sep =",")

# Plot coefficients
#===================
par(mar=c(6,6,3,6))
library(ggplot2)
coef <- as.data.frame(t(coef.matrix))
coef$vars <- rownames(coef)
ggplot(coef, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="log(svl)") +
  theme_classic()



#============================
# 8. PGSL: morpho ~ climatic
#============================
morpho <- read.csv("DATA/dt_morpho_55sp.csv")
morpho <- morpho[,-1]
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho <- morpho[morpho$sp!="Rana_luteiventris",]
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]
morpho <- morpho[morpho$sp!="Rana_aurora",]
morpho <- morpho[morpho$sp!="Rana_cascadae",]
rownames(morpho) <- morpho$sp

# Remove sp with NAs
morpho.red <- na.omit(morpho[,c(1:4,6:11)])  # 34 sp in total

# Calculate Residuals
mat.red <- NULL
for (i in 7:ncol(morpho.red)){
  res <- resid(lm(log(morpho.red[,i]) ~ log(morpho.red$SVL)))
  mat.red <- cbind(mat.red, res)
}

colnames(mat.red) <- c("res_HL", "res_HW", "res_FL", "res_TL")
rownames(mat.red) <- morpho.red$sp

# Merge datasets
mat.red <- cbind(morpho.red, mat.red) # 34 sp

# Clip phylogeny
library(castor)
phylo.red <- get_subtree_with_tips(phyloRana, mat.red$sp)
phylo.red <- phylo.red$subtree  # 34 sp
# Clip climatic dataset
clim.red <- clim3[mat.red$sp,]

# Same order as phylogeny
mat.red <- mat.red[phylo.red$tip.label,]
mat.red$sp==phylo.red$tip.label
rownames(clim.red)== phylo.red$tip.label
clim.red <- clim.red[phylo.red$tip.label,]

# covariance matrix GLOBAL
cov.mat.morpho <- vcv(phylo.red, model="Brownian") 

dt.clim <- rrpp.data.frame(clim.red)
dt.all <- rrpp.data.frame(clim = dt.clim, morpho = mat.red[,11:14], phy.cvc = cov.mat.morpho)
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=length(dt.all$morpho), ncol=length(dt.all[1:ncol(clim.red)]))

for (j in 1:ncol(dt.all$morpho)){
  for (i in 1:length(dt.all[1:ncol(clim.red)])){
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,i], y = dt.all$morpho[,j], phy.cvc = cov.mat.morpho)
    fitOLS <- lm.rrpp(y ~ x, data = dt.temp, print.progress = FALSE, iter=999)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
    ols <- cbind(anova(fitOLS)$table, colnames(clim.red)[i], colnames(dt.all$morpho)[j])
    gls <- cbind(anova(fitGLS)$table, colnames(clim.red)[i], colnames(dt.all$morpho)[j])
    
    anovas_OLS <- rbind(anovas_OLS, ols)
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim.red)[i], 
                                  colnames(dt.all$morpho)[j]))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[j,i] <- 0
    } else {
      coef.matrix[j,i] <- (fitGLS$LM$gls.coefficients[2])
    }
    
  }
}
colnames(coef.matrix) <- colnames(clim.red)
rownames(coef.matrix) <- colnames(dt.all$morpho)
coef.matrix <- round(coef.matrix, 3)

# SAVE
# write.table(anovas_OLS, "Results/PGLS/morpho/anovas_OLS.txt", sep = ",")
# write.table(anovas_GLS, "Results/PGLS/morpho/anovas_GLS.txt", sep = ",")
# write.table(coef_GLS, "Results/PGLS/morpho/coef_GLS.txt", sep = ",")
# write.table(coef.matrix, "Results/PGLS/morpho/coef_matrix.txt", sep =",")

# Plot coefficients
#===================
par(mar=c(6,6,3,6))

library(RColorBrewer)
library(pheatmap)


pheatmap(t(coef.matrix), cluster_cols = FALSE,
         cluster_rows = FALSE, border_color = "black", legend = F)

# Interactive heatmap
library(heatmaply)
heatmaply(t(coef.matrix), colors = colorRampPalette(brewer.pal(3, "RdBu"))(100),
          grid_gap = 1)

library("gplots")
col<- colorRampPalette(c("blue", "white", "red"))(100)
heatmap.2(t(coef.matrix), scale = "none", col = col, 
          trace = "none", density.info = "none", dendrogram= "none", 
          cexCol = 1, )

heatmap(t(coef.matrix))


# Separate by radiations
mat.red.EA <- mat.red[mat.red$phylo_reg=="EA",]
mat.red.AM <- mat.red[mat.red$phylo_reg=="AM",]

phylo.red.EA <- (get_subtree_with_tips(phyloRana, mat.red.EA$sp))$subtree  # 21 sp
phylo.red.AM <- (get_subtree_with_tips(phyloRana, mat.red.AM$sp))$subtree  # 13 sp

clim.red.EA <- clim.red[mat.red.EA$sp,]
clim.red.AM <- clim.red[mat.red.AM$sp,]


# 8.1 PGLS morpho EURASIA
#============================
cov.mat.morpho.EA <- vcv(phylo.red.EA, model="Brownian") 

dt.clim <- rrpp.data.frame(clim.red.EA)
dt.all <- rrpp.data.frame(clim = dt.clim, morpho = mat.red.EA[,11:14], phy.cvc = cov.mat.morpho.EA)
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=length(dt.all$morpho), ncol=length(dt.all[1:ncol(clim.red.EA)]))

for (j in 1:ncol(dt.all$morpho)){
  for (i in 1:length(dt.all[1:ncol(clim.red.EA)])){
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,i], y = dt.all$morpho[,j], phy.cvc = cov.mat.morpho.EA)
    fitOLS <- lm.rrpp(y ~ x, data = dt.temp, print.progress = FALSE, iter=999)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
    ols <- cbind(anova(fitOLS)$table, colnames(clim.red.EA)[i], colnames(dt.all$morpho)[j])
    gls <- cbind(anova(fitGLS)$table, colnames(clim.red.EA)[i], colnames(dt.all$morpho)[j])
    
    anovas_OLS <- rbind(anovas_OLS, ols)
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim.red.EA)[i], 
                                  colnames(dt.all$morpho)[j]))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[j,i] <- 0
    } else {
      coef.matrix[j,i] <- (fitGLS$LM$gls.coefficients[2])
    }
    
  }
}
colnames(coef.matrix) <- colnames(clim.red.EA)
rownames(coef.matrix) <- colnames(dt.all$morpho)
coef.matrix <- round(coef.matrix, 3)
# Any of the climatic variables have significant relationship with HL, HW, FL or TL in eurasian species


# 8.2 PGLS morpho AMERICA
#============================
cov.mat.morpho.AM <- vcv(phylo.red.AM, model="Brownian") 

dt.clim <- rrpp.data.frame(clim.red.AM)
dt.all <- rrpp.data.frame(clim = dt.clim, morpho = mat.red.AM[,11:14], phy.cvc = cov.mat.morpho.AM)
str(dt.all)

anovas_OLS  <- anovas_GLS <- NULL
coef_GLS <- NULL
coef.matrix <- matrix(nrow=length(dt.all$morpho), ncol=length(dt.all[1:ncol(clim.red.AM)]))

for (j in 1:ncol(dt.all$morpho)){
  for (i in 1:length(dt.all[1:ncol(clim.red.AM)])){
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,i], y = dt.all$morpho[,j], phy.cvc = cov.mat.morpho.AM)
    fitOLS <- lm.rrpp(y ~ x, data = dt.temp, print.progress = FALSE, iter=999)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.all$phy.cvc, print.progress = F)
    ols <- cbind(anova(fitOLS)$table, colnames(clim.red.AM)[i], colnames(dt.all$morpho)[j])
    gls <- cbind(anova(fitGLS)$table, colnames(clim.red.AM)[i], colnames(dt.all$morpho)[j])
    
    anovas_OLS <- rbind(anovas_OLS, ols)
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim.red.AM)[i], 
                                  colnames(dt.all$morpho)[j]))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[j,i] <- 0
    } else {
      coef.matrix[j,i] <- (fitGLS$LM$gls.coefficients[2])
    }
    
  }
}
colnames(coef.matrix) <- colnames(clim.red.AM)
rownames(coef.matrix) <- colnames(dt.all$morpho)
coef.matrix <- round(coef.matrix, 3)

# Plot coefficients
#===================
par(mar=c(8,8,8,8))
library('pheatmap')
pheatmap(t(coef.matrix), cluster_cols = FALSE,
         cluster_rows = FALSE, border_color = "black", legend = F)

# SAVE
# write.table(anovas_OLS, "Results/PGLS/morpho/anovas_OLS_AM.txt", sep = ",")
# write.table(anovas_GLS, "Results/PGLS/morpho/anovas_GLS_AM.txt", sep = ",")
# write.table(coef_GLS, "Results/PGLS/morpho/coef_GLS_AM.txt", sep = ",")
# write.table(coef.matrix, "Results/PGLS/morpho/coef_matrix_AM.txt", sep =",")

