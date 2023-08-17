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

phyloRana_68 <- read.tree("DATA/phyloRana_68sp.tre")
phyloEA <- read.tree("DATA/phyloEA_31sp.tre")
phyloAM <- read.tree("DATA/phyloAM_37sp.tre")


# 2. Import traits
#===================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
dt74 <- dt74[,-1]

dt68 <-read.csv("DATA/dt_68sp.csv")
dt68 <- dt68[,-1]
rownames(dt68) <- dt68$sp


# 3. Import climatic variables
#=============================
clim <- read.csv("DATA/clim_74sp.csv")
clim <- clim[,-c(1)]
rownames(clim) <- clim$X

# 4. Check order
#===============
phyloRana$tip.label == dt74$sp
dt74 <- dt74[phyloRana$tip.label,]

phyloRana_68$tip.label == dt68$sp
dt68 <- dt68[phyloRana_68$tip.label,]

phyloRana$tip.label == clim$X
clim <- clim[phyloRana$tip.label,]

clim68 <- clim[phyloRana_68$tip.label,]
phyloRana_68$tip.label == clim68$X


# 5. Scale climatic variables
#============================
# first remove SD
clim2 <- clim68[,-grep("sd", colnames(clim68))]
clim3 <- scale(clim2[,-1])


# 6. PGLS: lm.rrpp {SVL ~ all climatic variables}
#=================================================
# Build the Phylogenetic var-covariance matrix (correlation)
cov.mat.rana <- vcv(phyloRana_68, model="Brownian") 

dt.clim <- rrpp.data.frame(clim3)
dt.all <- rrpp.data.frame(clim = dt.clim, svl = dt68$log_svl, phy.cvc = cov.mat.rana)
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
coef.svl <- coef.matrix
par(mar=c(6,6,3,6))
library(ggplot2)
coef.svl <- as.data.frame(t(coef.svl))
coef.svl$vars <- rownames(coef.svl)
ggplot(coef.svl, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="log(svl)") +
  theme_classic()

# Cambiar los colores del plot
ggplot(coef.svl, aes(x=vars, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="coeffient") +
  theme_classic() + scale_fill_manual(values = c("#4F6367","#EEF5DB"),
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
coef.svl.ea <- coef.matrix
par(mar=c(6,6,3,6))
library(ggplot2)
coef.svl.ea <- as.data.frame(t(coef.svl.ea))
coef.svl.ea$vars <- rownames(coef.svl.ea)
ggplot(coef.svl.ea, aes(x=vars, y=svl)) +
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
coef.svl.am <- coef.matrix
par(mar=c(6,6,3,6))
library(ggplot2)
coef.svl.am <- as.data.frame(t(coef.svl.am))
coef.svl.am$vars <- rownames(coef.svl.am)
ggplot(coef.svl.am, aes(x=vars, y=svl)) +
  geom_bar(stat="identity", aes(fill = ifelse(svl <0, "red", "blue"))) + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="log(svl)") +
  theme_classic()

# PLOT all together
#===================
p1 <- ggplot(coef.svl, aes(x=vars, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="coefficient") +
  theme_classic() + 
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("Global") +
  geom_hline(yintercept = 0, color="gray50")

p2 <- ggplot(coef.svl.ea, aes(x=vars, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="coefficient") +
  theme_classic() + 
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("Eurasia") +
  geom_hline(yintercept = 0, color="gray50")

p3 <- ggplot(coef.svl.am, aes(x=vars, y=svl, fill=svl > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + coord_flip()+labs(x= "", y="coefficient") +
  theme_classic() + 
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("America") +
  geom_hline(yintercept = 0, color="gray50")


library(patchwork)
p1 + p2 + p3



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


# MEGALOOP to not miss species and clip dataset and phylogeny by each trait
library(castor)
morpho.res <- morpho[, c(1,13:16)]
mat <- NULL

for (i in 2:ncol(morpho.res)){
  morpho.x <- morpho.res[complete.cases(morpho.res[,i]), c(1,i)]      # Remove sp from morpho dataset
  clim.x <- clim3[morpho.x$sp,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.x$sp))$subtree  # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- clim.x[phylo.x$tip.label,]
  
  cov.temp <- vcv(phylo.x, model="Brownian") 
  
  dt.clim <- rrpp.data.frame(clim.x)
  dt.all <- rrpp.data.frame(clim = dt.clim, morpho = morpho.x[,2], phy.cvc = cov.temp)
  
  anovas_GLS <- NULL
  coef_GLS <- NULL
  coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(clim.x)]))
  
  for (j in 1:length(dt.all[1:ncol(clim.x)])){   # The normal loop for y ~ each climatic variable
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,j], y = dt.all$morpho, phy.cvc = cov.temp)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.temp$phy.cvc, print.progress = F)
    gls <- cbind(anova(fitGLS)$table, colnames(clim.x)[j], colnames(morpho.x)[2])
    
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim.x)[j], 
                                  colnames(morpho.x)[2]))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[1,j] <- 0
    } else {
      coef.matrix[1,j] <- (fitGLS$LM$gls.coefficients[2])
    }
    
    colnames(coef.matrix) <- colnames(clim.x)
    rownames(coef.matrix) <- colnames(morpho.x)[2]
    coef.matrix <- round(coef.matrix, 3)
  }
  mat <- rbind(mat, coef.matrix)

}

# Plot coefficients
library("gplots")
library(pheatmap)
library(RColorBrewer)
col<- colorRampPalette(c("tomato", "white", "steelblue3"))(100)
heatmap.2(t(mat), scale = "none", col = col, 
          trace = "none", density.info = "none", dendrogram= "none", 
          cexCol = 1, Colv=FALSE, Rowv = FALSE)
pheatmap(t(mat), cluster_cols = FALSE,
         cluster_rows = FALSE, border_color = "black", legend = F,
         color = brewer.pal(n = 11, name = "RdBu"))




# Separate by radiations
#=========================

##########
# EURASIA
##########
morpho.ea <- morpho[morpho$phylo_reg=="EA",]
morpho.ea <- morpho.ea[, c(1,13:16)]
mat.ea <- NULL
all.gls <- NULL
all.coefs <- NULL

for (i in 2:ncol(morpho.ea)){
  morpho.x <- morpho.ea[complete.cases(morpho.ea[,i]), c(1,i)]      # Remove sp from morpho dataset
  clim.x <- clim3[morpho.ea$sp,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.ea$sp))$subtree  # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- clim.x[phylo.x$tip.label,]
  
  cov.temp <- vcv(phylo.x, model="Brownian") 
  
  dt.clim <- rrpp.data.frame(clim.x)
  dt.all <- rrpp.data.frame(clim = dt.clim, morpho = morpho.x[,2], phy.cvc = cov.temp)
  
  anovas_GLS <- NULL
  coef_GLS <- NULL
  coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(clim.x)]))
  
  for (j in 1:length(dt.all[1:ncol(clim.x)])){   # The normal loop for y ~ each climatic variable
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,j], y = dt.all$morpho, phy.cvc = cov.temp)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.temp$phy.cvc, print.progress = F)
    gls <- cbind(anova(fitGLS)$table, colnames(clim.x)[j], colnames(morpho.x)[2])
    
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim.x)[j], 
                                  colnames(morpho.x)[2]))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[1,j] <- 0
    } else {
      coef.matrix[1,j] <- (fitGLS$LM$gls.coefficients[2])
    }
    
    colnames(coef.matrix) <- colnames(clim.x)
    rownames(coef.matrix) <- colnames(morpho.x)[2]
    coef.matrix <- round(coef.matrix, 3)
  }
  mat.ea <- rbind(mat.ea, coef.matrix)
  all.gls <- rbind(all.gls, anovas_GLS)
  all.coefs <- rbind(all.coefs, coef_GLS)
  
}  

# SAVE
# write.table(all.gls, "Results/PGLS/morpho/anovas_GLS_EA.txt", sep = ",")
# write.table(all.coefs, "Results/PGLS/morpho/coef_GLS_EA.txt", sep = ",")
# write.table(mat.ea, "Results/PGLS/morpho/coef_matrix_EA.txt", sep =",")


# Plot coefficients
col<- colorRampPalette(c("tomato", "white", "steelblue3"))(100)
heatmap.2(t(mat.ea), scale = "none", col = col, 
          trace = "none", density.info = "none", dendrogram= "none", 
          cexCol = 1, Colv=FALSE, Rowv = FALSE)
pheatmap(t(mat.ea), cluster_cols = FALSE,
         cluster_rows = FALSE, border_color = "black", legend = F,
         color = brewer.pal(n = 11, name = "RdBu"))



##########
# AMERICA
##########
morpho.am <- morpho[morpho$phylo_reg=="AM",]
morpho.am <- morpho.am[, c(1,13:16)]
mat.am <- NULL
all.gls.am <- NULL
all.coefs.am <- NULL

for (i in 2:ncol(morpho.am)){
  morpho.x <- morpho.am[complete.cases(morpho.am[,i]), c(1,i)]      # Remove sp from morpho dataset
  clim.x <- clim3[morpho.am$sp,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.am$sp))$subtree  # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- clim.x[phylo.x$tip.label,]
  
  cov.temp <- vcv(phylo.x, model="Brownian") 
  
  dt.clim <- rrpp.data.frame(clim.x)
  dt.all <- rrpp.data.frame(clim = dt.clim, morpho = morpho.x[,2], phy.cvc = cov.temp)
  
  anovas_GLS <- NULL
  coef_GLS <- NULL
  coef.matrix <- matrix(nrow=1, ncol=length(dt.all[1:ncol(clim.x)]))
  
  for (j in 1:length(dt.all[1:ncol(clim.x)])){   # The normal loop for y ~ each climatic variable
    dt.temp <- rrpp.data.frame(x = dt.all[[1]][,j], y = dt.all$morpho, phy.cvc = cov.temp)
    fitGLS <- lm.rrpp(y ~ x, data = dt.temp, Cov = dt.temp$phy.cvc, print.progress = F)
    gls <- cbind(anova(fitGLS)$table, colnames(clim.x)[j], colnames(morpho.x)[2])
    
    anovas_GLS <- rbind(anovas_GLS, gls)
    
    coef_GLS <- rbind(coef_GLS, c(fitGLS$LM$gls.coefficients, colnames(clim.x)[j], 
                                  colnames(morpho.x)[2]))
    
    if (anova(fitGLS)$table[1,7] > 0.05) {
      coef.matrix[1,j] <- 0
    } else {
      coef.matrix[1,j] <- (fitGLS$LM$gls.coefficients[2])
    }
    
    colnames(coef.matrix) <- colnames(clim.x)
    rownames(coef.matrix) <- colnames(morpho.x)[2]
    coef.matrix <- round(coef.matrix, 3)
  }
  mat.am <- rbind(mat.am, coef.matrix)
  all.gls.am <- rbind(all.gls.am, anovas_GLS)
  all.coefs.am <- rbind(all.coefs.am, coef_GLS)
  
}  

# SAVE
# write.table(all.gls.am, "Results/PGLS/morpho/anovas_GLS_AM.txt", sep = ",")
# write.table(all.coefs.am, "Results/PGLS/morpho/coef_GLS_AM.txt", sep = ",")
# write.table(mat.am, "Results/PGLS/morpho/coef_matrix_AM.txt", sep =",")


# Plot coefficients
col<- colorRampPalette(c("tomato", "white", "steelblue3"))(100)
heatmap.2(t(mat.am), scale = "none", col = col, 
          trace = "none", density.info = "none", dendrogram= "none", 
          cexCol = 1, Colv=FALSE, Rowv = FALSE, main="American")
pheatmap(t(mat.am), cluster_cols = FALSE,
         cluster_rows = FALSE, border_color = "black", legend = F,
         color = brewer.pal(n = 11, name = "RdBu"))




# PLOT all together
#===================
library(ComplexHeatmap)
library(circlize)
library(tiff)

round(seq(-0.3, 0.3, length.out = 9),3)
paleta <- colorRamp2(c(-0.300, -0.225, -0.150, -0.075,  0.000,  0.075,  0.150,  0.225,  0.300), 
                     c("#67001F", "#B2182B", "#D6604D",  "#F4A582", "#FFFFFF", "#92C5DE",
                       "#4393C3",  "#2166AC", "#053061"))

ht_opt(heatmap_column_names_gp = gpar(fontsize = 9),
       heatmap_row_names_gp = gpar(fontsize = 9),
       heatmap_column_title_gp = gpar(fontsize = 16, fontface="italic"),
       legend_grid_height = unit(10, "cm"),
       legend_labels_gp = gpar(fonsize=8, cex = 1),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = FALSE)

ht1 <- Heatmap(t(mat), col = paleta, cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               row_names_side = "left", border = T, heatmap_height = unit(15, "cm"), 
               column_title = "Global", name=" ", rect_gp = gpar(col="#F4F4F4")) 

ht2 <- Heatmap(t(mat.ea), col = paleta, cluster_rows = FALSE, 
               cluster_columns = FALSE, show_heatmap_legend = FALSE,
               row_names_side = "left", border = T, heatmap_height = unit(15, "cm"), 
               column_title = "Eurasia", name=" ", rect_gp = gpar(col="#F4F4F4")) 

ht3 <- Heatmap(t(mat.am), col = paleta, cluster_rows = FALSE, 
               cluster_columns = FALSE, show_heatmap_legend = FALSE,
               row_names_side = "left", border = T, heatmap_height = unit(15, "cm"), 
               column_title = "America", name=" ", rect_gp = gpar(col="#F4F4F4"))

ht_list = ht1 + ht2 + ht3

par(mar=c(1,1,1,1), oma=c(1,1,1,1), cex.axis = 0.5, font.axis = 2)
draw(ht_list, ht_gap = unit(0.5, "cm"))


