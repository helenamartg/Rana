#################################
# Phylogenetic linear models
# with BM and OU: 
# MORPHOLOGICAL VARIABLES
#################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

library(phylolm)
library(phytools)
library(castor)


# 1. Import traits: morphological traits 
#===========================================
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

# 2. Import phylogeny
#======================
phyloEA <- read.tree("Data/phyloEA_26sp.tre")
phyloAM <- read.tree("Data/phyloAM_33sp.tre")

phyloEA$root.edge <- NULL
phyloAM$root.edge <- NULL


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
rownames(climEA)==phyloEA$tip.label
rownames(climAM)==phyloAM$tip.label

climEA <- climEA[phyloEA$tip.label,]
climAM <- climAM[phyloAM$tip.label,]



# 5. Allen's rule:  TMINavg + SRADavg
#=================

# 5.1. Eurasia: 
###############
morpho.ea <- morpho[morpho$phylo_reg=="EA",]
morpho.ea <- morpho.ea[, c(1,12:15)]


results_morpho_ea <- data.frame()
mat.bm.ea <- data.frame()
mat.ou.ea <- data.frame()

for(i in 2:ncol(morpho.ea)) {
  morpho.x <- morpho.ea[complete.cases(morpho.ea[,i]), c(1,i)]          # Remove sp from morpho dataset
  clim.x <- climEA[morpho.ea$Species,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloEA, morpho.x$Species))$subtree      # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- as.data.frame(clim.x[phylo.x$tip.label,])
  
  trait <- setNames(morpho.x[,2], morpho.x$Species)
  
  # Set the root edge to 0 (or remove it entirely)
  if(!is.null(phylo.x$root.edge)) {
    phylo.x$root.edge <- NULL  # Remove root edge if it exists
  }

  bm_model_eurasia <- phylolm(trait ~ TMINavg+ SRADavg, data=clim.x,
                                phy = phylo.x, model = "BM")
  ou_model_eurasia <- phylolm(trait ~ TMINavg+ SRADavg, data=clim.x,
                                phy = phylo.x, model = "OUfixedRoot")
    
  bm_summary <- summary(bm_model_eurasia)
  ou_summary <- summary(ou_model_eurasia)
    
  bm <- as.data.frame(cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                              Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik))
  bm$trait <- as.character(colnames(morpho.x)[2])
  rownames(bm)[2] <- c("TMINavg+ SRADavg")
    
  ou <- as.data.frame(cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                              Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik))
  ou$trait <- as.character(colnames(morpho.x)[2])
  rownames(ou)[2] <- c("TMINavg+ SRADavg")
    
  mat.bm.ea <- rbind(mat.bm.ea, bm)
  mat.ou.ea <- rbind(mat.ou.ea, ou)
    
  AICc_bm <- AIC(bm_model_eurasia, corrected = TRUE)
  AICc_ou <- AIC(ou_model_eurasia, corrected = TRUE)
    
  x <- cbind(trait = colnames(morpho.x)[2], 
               clim_var = c("TMINavg+ SRADavg"), 
               AICc_BM = as.numeric(AICc_bm), AICc_OU = as.numeric(AICc_ou),
               Radiation = "Eurasia", Alpha = ou_summary$optpar, Lowest_AICc=NA)
    
  results_morpho_ea <- rbind(results_morpho_ea, x)
    
}

results_morpho_ea <- as.data.frame(results_morpho_ea)
results_morpho_ea[,3] <- as.numeric(results_morpho_ea[,3])
results_morpho_ea[,4] <- as.numeric(results_morpho_ea[,4])
results_morpho_ea[,6] <- as.numeric(results_morpho_ea[,6])



# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results_morpho_ea)){
  if (results_morpho_ea[i,]$AICc_BM < results_morpho_ea[i,]$AICc_OU) {
    results_morpho_ea[i,]$Lowest_AICc <- "BM model"
  } else {
    results_morpho_ea[i,]$Lowest_AICc <- "OU model"
  }
}


# Calculate delta AICc
library(dplyr)
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
  if (results_morpho_ea[i,11]==TRUE){
    results_morpho_ea[i,12] <- results_morpho_ea[i,7]
  } else {
    results_morpho_ea[i,12] <- "BM"
  }
}

print(results_morpho_ea)






# 5.2. America: 
###############
morpho.am <- morpho[morpho$phylo_reg=="AM",]
morpho.am <- morpho.am[, c(1,12:15)]


results_morpho_am <- data.frame()
mat.bm.am <- data.frame()
mat.ou.am <- data.frame()

for(i in 2:ncol(morpho.am)) {
  morpho.x <- morpho.am[complete.cases(morpho.am[,i]), c(1,i)]          # Remove sp from morpho dataset
  clim.x <- climAM[morpho.am$Species,]                                       # Clip climatic dataset
  phylo.x <- (get_subtree_with_tips(phyloAM, morpho.x$Species))$subtree      # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  clim.x <- as.data.frame(clim.x[phylo.x$tip.label,])
  
  trait <- setNames(morpho.x[,2], morpho.x$Species)
  
  # Set the root edge to 0 (or remove it entirely)
  if(!is.null(phylo.x$root.edge)) {
    phylo.x$root.edge <- NULL  # Remove root edge if it exists
  }
  
  bm_model_america <- phylolm(trait ~ TMINavg + SRADavg, data=clim.x,
                              phy = phylo.x, model = "BM")
  ou_model_america <- phylolm(trait ~ TMINavg + SRADavg, data=clim.x,
                              phy = phylo.x, model = "OUfixedRoot")
  
  bm_summary <- summary(bm_model_america)
  ou_summary <- summary(ou_model_america)
  
  bm <- as.data.frame(cbind(bm_summary$coefficients, Rsq = bm_summary$r.squared, 
                            Rsq_adj = bm_summary$adj.r.squared, logLik = bm_summary$logLik))
  bm$trait <- as.character(colnames(morpho.x)[2])
  rownames(bm)[2] <- c("TMINavg+ SRADavg")
  
  ou <- as.data.frame(cbind(ou_summary$coefficients, Rsq = ou_summary$r.squared, 
                            Rsq_adj = ou_summary$adj.r.squared, logLik = ou_summary$logLik))
  ou$trait <- as.character(colnames(morpho.x)[2])
  rownames(ou)[2] <- c("TMINavg+ SRADavg")
  
  mat.bm.am <- rbind(mat.bm.am, bm)
  mat.ou.am <- rbind(mat.ou.am, ou)
  
  AICc_bm <- AIC(bm_model_america, corrected = TRUE)
  AICc_ou <- AIC(ou_model_america, corrected = TRUE)
  
  x <- cbind(trait = colnames(morpho.x)[2], 
             clim_var = c("TMINavg+ SRADavg"), 
             AICc_BM = as.numeric(AICc_bm), AICc_OU = as.numeric(AICc_ou),
             Radiation = "America", Alpha = ou_summary$optpar, Lowest_AICc=NA)
  
  results_morpho_am <- rbind(results_morpho_am, x)
  
}

results_morpho_am <- as.data.frame(results_morpho_am)
results_morpho_am[,3] <- as.numeric(results_morpho_am[,3])
results_morpho_am[,4] <- as.numeric(results_morpho_am[,4])
results_morpho_am[,6] <- as.numeric(results_morpho_am[,6])



# Compare AICcs of BM and OU models of each climatic variable
for (i in 1:nrow(results_morpho_am)){
  if (results_morpho_am[i,]$AICc_BM < results_morpho_am[i,]$AICc_OU) {
    results_morpho_am[i,]$Lowest_AICc <- "BM model"
  } else {
    results_morpho_am[i,]$Lowest_AICc <- "OU model"
  }
}


# Calculate delta AICc
library(dplyr)
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
  if (results_morpho_am[i,11]==TRUE){
    results_morpho_am[i,12] <- results_morpho_am[i,7]
  } else {
    results_morpho_am[i,12] <- "BM"
  }
}

print(results_morpho_am)



# 6. MERGE RADIATIONS
#======================

results <- rbind(results_morpho_ea, results_morpho_am)

mat.bm.am$model <- "BM"
mat.ou.am$model <- "OU"
mat.am <- rbind(mat.bm.am, mat.ou.am)

mat.bm.ea$model <- "BM"
mat.ou.ea$model <- "OU"
mat.ea <- rbind(mat.bm.ea, mat.ou.ea)


# SAVE
# write.csv(results, "Results/phylolm/morpho/Results_both_radiations.csv")
# write.csv(mat.ea, "Results/phylolm/morpho/coefmat_Eurasia.csv")
# write.csv(mat.am, "Results/phylolm/morpho/coefmat_America.csv")





# PLOT morpho ~ climate 
##########################

# EURASIA: 
mat.bm.ea
eurasia <- mat.bm.ea[,c(1,4,8)]
eurasia <- cbind(rownames(eurasia), eurasia)
names(eurasia) <- c("clim", "Estimate","p-value", "trait")

coef.mat <- matrix(nrow=nrow(eurasia), ncol = 3)
for( i in 1:nrow(eurasia)) {
  coef.mat[i,1] <- eurasia[i,1]
  coef.mat[i,2] <- eurasia[i,4]
  if (eurasia$`p-value`[i] > 0.05) {
    coef.mat[i,3] <- 0
  } else {
    coef.mat[i,3] <- round((eurasia$Estimate[i]),3)
  }
  
}

coef.mat <- as.data.frame(coef.mat)
colnames(coef.mat) <- c("clim.var", "trait", "estimate")
rownames(coef.mat) <- coef.mat$clim.var
coef.mat$estimate <- as.numeric(coef.mat$estimate)

HL <- (coef.mat[coef.mat$trait=="phylores_HL",])[,-1]
HW <- (coef.mat[coef.mat$trait=="phylores_HW",])[,-1]
FL <- (coef.mat[coef.mat$trait=="phylores_FL",])[,-1]
TL <- (coef.mat[coef.mat$trait=="phylores_TL",])[,-1]

coef.mat2 <- rbind(t(HL[,-1]), t(HW[,-1]), t(FL[,-1]), t(TL[,-1])) 
rownames(coef.mat2) <- unique(coef.mat$trait)
coef.mat2 <- coef.mat2[,-1]
colnames(coef.mat2) <- c("TMINavg", "SRADavg")



# AMERICA: 
# HL: all BM
# HW: all OU
# FL: all BM
# TL: all BM


hl <-  mat.bm.am[mat.bm.am$trait=="phylores_HL",]
hw <- mat.ou.am[mat.ou.am$trait=="phylores_HW",]

fl <- mat.bm.am[mat.bm.am$trait=="phylores_FL",]
tl <- mat.bm.am[mat.bm.am$trait=="phylores_TL",]


mat.bm.def <- rbind(hl, hw, fl, tl)
america <- mat.bm.def[,c(1,4,8)]
america <- cbind(rownames(america), america)
names(america) <- c("clim", "Estimate","p-value", "trait")

coef.mat <- matrix(nrow=nrow(america), ncol = 3)
for( i in 1:nrow(america)) {
  coef.mat[i,1] <- america[i,1]
  coef.mat[i,2] <- america[i,4]
  if (america$`p-value`[i] > 0.05) {
    coef.mat[i,3] <- 0
  } else {
    coef.mat[i,3] <- round((america$Estimate[i]),3)
  }
  
}

coef.mat <- as.data.frame(coef.mat)
colnames(coef.mat) <- c("clim.var", "trait", "estimate")
rownames(coef.mat) <- coef.mat$clim.var
coef.mat$estimate <- as.numeric(coef.mat$estimate)

HL <- (coef.mat[coef.mat$trait=="phylores_HL",])[,-1]
HW <- (coef.mat[coef.mat$trait=="phylores_HW",])[,-1]
FL <- (coef.mat[coef.mat$trait=="phylores_FL",])[,-1]
TL <- (coef.mat[coef.mat$trait=="phylores_TL",])[,-1]

coef.mat3 <- rbind(t(HL[,-1]), t(HW[,-1]), t(FL[,-1]), t(TL[,-1])) 
rownames(coef.mat3) <- unique(coef.mat$trait)
coef.mat3 <- coef.mat3[,-1]
colnames(coef.mat3) <- c("TMINavg", "SRADavg")



# Plot coefficients
library(gplots)
library(pheatmap)
library(ComplexHeatmap)

# PLOT all together
#===================
library(circlize)
library(tiff)

round(seq(-0.3, 0.3, length.out = 9),3)
paleta <- colorRamp2(c(-0.300, -0.225, -0.150, -0.075,  0.000,  0.075,  0.150,  0.225,  0.300), 
                     c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#FFFFFF",
                       "#F4A582", "#D6604D", "#B2182B","#67001F"))


   
   


library("ComplexHeatmap")
ht_opt(heatmap_column_names_gp = gpar(fontsize = 9),
       heatmap_row_names_gp = gpar(fontsize = 9),
       heatmap_column_title_gp = gpar(fontsize = 16, fontface="italic"),
       legend_grid_height = unit(10, "cm"),
       legend_labels_gp = gpar(fonsize=8, cex = 1),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = FALSE)


ht2 <- Heatmap(t(coef.mat2), col = paleta, cluster_rows = FALSE, 
               cluster_columns = FALSE, show_heatmap_legend = FALSE,
               row_names_side = "left", border = T, heatmap_height = unit(15, "cm"), 
               column_title = "Eurasia", name=" ", rect_gp = gpar(col="#F4F4F4")) 

ht3 <- Heatmap(t(coef.mat3), col = paleta, cluster_rows = FALSE, 
               cluster_columns = FALSE, show_heatmap_legend = FALSE,
               row_names_side = "left", border = T, heatmap_height = unit(15, "cm"), 
               column_title = "America", name=" ", rect_gp = gpar(col="#F4F4F4"))

ht_list = ht2 + ht3

par(mar=c(1,1,1,1), oma=c(1,1,1,1), cex.axis = 0.5, font.axis = 2)
draw(ht_list, ht_gap = unit(0.5, "cm"))






# Barplot coefficients
#######################
HL.ea <- as.data.frame(cbind(as.numeric(coef.mat2[1,]), c("TMINavg", "SRADavg")))
colnames(HL.ea) <- c("Estimate", "Variable")
HL.ea$Estimate <- as.numeric(as.character(HL.ea$Estimate))

p <- ggplot(HL.ea, aes(x=Variable, y=Estimate, fill=Estimate > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + labs(x= "", y="coefficient") +
  theme_classic() + 
  coord_flip() +
  # scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("Eurasia") +
  geom_hline(yintercept = 0, color="gray50")+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=15))


# America
HL.am <- as.data.frame(cbind(as.numeric(coef.mat3[1,]), c("TMINavg", "SRADavg")))
colnames(HL.am) <- c("Estimate", "Variable")
HL.am$Estimate <- as.numeric(as.character(HL.am$Estimate))

q <- ggplot(HL.am, aes(x=Variable, y=Estimate, fill=Estimate > 0)) +
  geom_bar(stat="identity") + 
  theme(legend.position = 'none') + labs(x= "", y="coefficient") +
  theme_classic() + 
  coord_flip() +
  # scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("#DD535E","steelblue3"),
                    labels=c("", ""))+
  labs(fill='svl') + ylim(-0.5,0.5) + guides(fill = FALSE) + ggtitle("America") +
  geom_hline(yintercept = 0, color="gray50")+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=15))


# Arrange in two columns
library(gridExtra)
grid.arrange(p, q, ncol=2)


