#############################
#   Cleaning morpho Ranidae
#############################
rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4publication/2nd_Revision/Reanalisis")

# 1. Import data
#=================
svl <- read.csv("Data/body_size_Rana.csv", sep=";")
rownames(svl) <- svl$Species

mor <- read.csv("Data/Morpho_rana.csv", sep=";")
mor <- mor[,-2]
rownames(mor) <- mor$Species

for (i in 1:nrow(svl)){
  if (svl[i,2]=="Lithobates"){
    svl[i,1] <- gsub("^Lithobates_", "Rana_", svl[i,1])
  }
}

# remove species for which we don't have SVL max
svl$SVL.max..def.
svl <- svl[-c(66:73),]


# 2. Combine datasets
#=====================
# Combinar los datasets usando la columna común (Species)
# dt <- merge(svl, mor, by = "Species", all.x = TRUE)

dt <- svl

# 3. Clean de dataset
#========================
names(dt)
dt <- dt[,-c(7:10)]
unique(dt$ECOTYPE)
# dt_red <- dt[,-c(6:26)]

# Remove Rana kunyensis (= R. coreana)
dt <- dt[dt$Species!="Rana_kunyuensis",]

# Remove Rana tlaloci (no climatic data)
dt <- dt[dt$Species!="Rana_tlaloci",]

# Remove Rana psilonota (no SVL data)
dt <- dt[dt$Species!="Rana_psilonota",]

dt65 <- dt

# Add the real region
rownames(dt65) <- dt65$Species

dt65$real_reg <- NA
for (i in 1:nrow(dt65)){
  if (dt65[i,4]=="AS" | dt65[i,4]=="EU"){
    dt65[i,7] <- "EA"
  } else {
    dt65[i,7] <- "AM"
  }
}

# Reorder
names(dt65)
dt65 <- dt65[, c(1:4, 7, 5,6)]


# Log-transformed body size 
dt65$log_svl <- log(dt65$SVL.max..def.)


# SAVE the reduced dataset with 65 sp 
# write.csv(dt65, "Data/dt_65sp.csv")

tr74 <- read.tree("DATA/phyloRana_74sp.tre")
tr65 <- (get_subtree_with_tips(tr74, dt65$Species))$subtree 

# Save tree
# write.tree(tr65, "DATA/phyloRana_65sp.tre")



# 4. Obtain a dataset of morphological variables
#=================================================
morpho <- mor
rownames(morpho)

# Missing values
apply(morpho, 2, function(x){sum(is.na(x))})

# Leave only morphological variables : HL, HW, FL, TL
def_morph <- morpho[,c(1,2, 4, 5, 18, 20)]
def_morph <- as.data.frame(def_morph)

# Remove species for which morphological data is not enough
def_morph$Species
def_morph <- def_morph[def_morph$Species!="Rana_neovolcanica",]
def_morph <- def_morph[def_morph$Species!="Rana_pustulosa",]
def_morph <- def_morph[def_morph$Species!="Rana_spectabilis",]
def_morph <- def_morph[def_morph$Species!="Rana_zweifeli",]
def_morph <- def_morph[def_morph$Species!="Rana_asiatica",]
def_morph <- def_morph[def_morph$Species!="Rana_zhenhaiensis",]
def_morph <- def_morph[def_morph$Species!="Rana_hanluica",]
def_morph <- def_morph[def_morph$Species!="Rana_chaochiaoensis",]
def_morph <- def_morph[def_morph$Species!="Rana_weiningensis",]
def_morph <- def_morph[def_morph$Species!="Rana_tlaloci",]

apply(def_morph, 2, function(x){sum(is.na(x))}) # HW for all the 47 sp

# Reorder
def_morph$Species <- rownames(def_morph)
names(def_morph)

def_morph2 <- merge(dt, def_morph, by = "Species")
names(def_morph2)
def_morph2 <- def_morph2[,-7]

colnames(def_morph2) <- c("Species", "genus", "phylo_reg", "state", "ecotype",
                    "SVLmax", "HL", "HW","FL", "TL")

# log-transformed
dt45 <- as.data.frame(def_morph2)
dt45$log_svl <- log(as.numeric(dt45$SVL))



################################
# phylogenetic size-correction
################################
library(ape)
library(phytools)
library(castor)

rownames(dt45) <- dt45$Species ## 45
tr55 <- read.tree("DATA/phyloRana_55sp.tre")

tr45 <- (get_subtree_with_tips(tr55, dt45$Species))$subtree 

# Save tree
write.tree(tr45, "DATA/phyloRana_45sp.tre")

# Check the order
rownames(dt45) == tr45$tip.label
dt45 <- dt45[tr45$tip.label,]


svl <- setNames(as.numeric(dt45$log_svl), dt45$Species)
hw <- setNames(log(as.numeric(dt45$HW)), dt45$Species)
hl <- setNames(log(as.numeric(dt45$HL)), dt45$Species)
fl <- setNames(log(as.numeric(dt45$FL)), dt45$Species)
tl <- setNames(log(as.numeric(dt45$TL)), dt45$Species)

class(dt45)

# Remove NAs
hl <- na.omit(dt45[,c(1,7,11)])
hw <-na.omit(dt45[,c(1,8,11)])
fl <-na.omit(dt45[,c(1,9,11)])
tl <- na.omit(dt45[,c(1,10,11)])

# Calculate phylogenetic residuals
list <- list(hl, hw, fl, tl)
phylores.list <- list(NULL, NULL, NULL, NULL)

for (i in 1:length(list)){
  svl <- setNames(as.numeric(list[[i]]$log_svl), list[[i]]$Species)
  morph <- setNames(log(as.numeric(list[[i]][,2])), list[[i]]$Species)
  tr <- (get_subtree_with_tips(tr45, list[[i]]$Species))$subtree  # clip phylogeny
  
  p <- phyl.resid(tr, svl, morph)$resid
  colnames(p) <-  paste("phylores", colnames(list[[i]][2]), sep = "_")
  
  phylores.list[[i]] <- p

}


# put all in a matrix
phylores.dt <- as.data.frame(matrix(nrow=length(dt45$Species)))
rownames(phylores.dt) <- dt45$Species
phylores.dt[,1] <- as.factor(rownames(phylores.dt))



# HL
phylores.dt <- cbind(phylores.dt, NA)

for (j in 1:nrow(phylores.list[[1]])){
  for (k in 1:nrow(phylores.dt)){
    if (rownames(phylores.list[[1]])[j] == phylores.dt[k,1]){
      x <- phylores.list[[1]][j]
      phylores.dt[k,2] <- as.numeric(x)
    } 
  }
}


# HW
phylores.dt <- cbind(phylores.dt, NA)

for (j in 1:nrow(phylores.list[[2]])){
  for (k in 1:nrow(phylores.dt)){
    if (rownames(phylores.list[[2]])[j] == phylores.dt[k,1]){
      x <- phylores.list[[2]][j]
      phylores.dt[k,3] <- as.numeric(x)
    } 
  }
}


# FL
phylores.dt <- cbind(phylores.dt, NA)

for (j in 1:nrow(phylores.list[[3]])){
  for (k in 1:nrow(phylores.dt)){
    if (rownames(phylores.list[[3]])[j] == phylores.dt[k,1]){
      x <- phylores.list[[3]][j]
      phylores.dt[k,4] <- as.numeric(x)
    } 
  }
}


# TL
phylores.dt <- cbind(phylores.dt, NA)

for (j in 1:nrow(phylores.list[[4]])){
  for (k in 1:nrow(phylores.dt)){
    if (rownames(phylores.list[[4]])[j] == phylores.dt[k,1]){
      x <- phylores.list[[4]][j]
      phylores.dt[k,5] <- as.numeric(x)
    } 
  }
}


dim(phylores.dt)
phylores.dt <- as.data.frame(phylores.dt)
colnames(phylores.dt) <- c("Species", paste("phylores", colnames(dt45[,7:10]), sep = "_"))


# merge with original dataset
dt45$Species == phylores.dt$Species
dt45 <- cbind(dt45, phylores.dt[,-1])

class(dt45)


# Calculate correlation between old residuals and phylo.residuals
cor(as.numeric(dt45$res_HL), as.numeric(dt45$phylores_HL), method= "pearson")
cor(as.numeric(dt45$res_HW), as.numeric(dt45$phylores_HW), use= "pairwise.complete.obs")
cor(dt45$res_FL, as.numeric(dt45$phylores_FL), use= "pairwise.complete.obs")
cor(dt45$res_TL, as.numeric(dt45$phylores_TL), use= "pairwise.complete.obs")


# Because correlations are very high (> 0.98 with the exception of HL 0.88), 
# our results shouldn't change if we perform new analyses with this phylogenetic correcction.

# SAVE
write.csv(dt45, "Data/morpho_phylores_45sp.csv")


# remove 6 species from Eurasian radiation distributed in NAM
dt65 <- dt65[dt65$Species!="Rana_luteiventris",]
dt65 <- dt65[dt65$Species!="Rana_boylii",]
dt65 <- dt65[dt65$Species!="Rana_sierrae",]
dt65 <- dt65[dt65$Species!="Rana_muscosa",]
dt65 <- dt65[dt65$Species!="Rana_aurora",]
dt65 <- dt65[dt65$Species!="Rana_cascadae",]
dt59 <- dt65

tr65 <- read.tree("Data/phyloRana_65sp.tre")
tr59 <- (get_subtree_with_tips(tr65, dt59$Species))$subtree 

# # SAVE
# write.csv(dt59, "Data/dt_59sp.csv")
# write.tree(tr59, "Data/phyloRana_59sp.tre")




# Separate radiations
#======================
tr65  # all species with SVL max
tr59  # Without 6 eurasian species
plot(tr59, cex=0.6)
nodelabels(cex=0.6)


phyloAM <- extract.clade(tr59, 62)
plot(phyloAM)
phyloEA <- extract.clade(tr59, 61)
plot(phyloEA)

plot(phyloAM)
is.ultrametric(phyloAM)

plot(phyloEA)
is.ultrametric(phyloEA)

table(dt59$real_reg)
# AM EA 
# 33 26 


# SAVE
# write.tree(phyloEA, "Data/phyloEA_26sp.tre")
# write.tree(phyloAM, "Data/phyloAM_33sp.tre")









#This is extra because of the reviewer's comment
##############################################################################
# Phylogenetic Residuals: Body Size Correction
# Intercept-only model

tr74 <- read.tree("DATA/phyloRana_74sp.tre")
dt74 <- read.csv("DATA/dt_74sp.csv")
rownames(dt74) <- dt74$sp

dt74$sp == tr74$tip.label
dt74 <- dt74[tr74$tip.label,]


library(caper)
# Prepara los datos para el análisis PGLS
comparative_data <- comparative.data(phy = tr74, data = dt74, names.col = "sp", vcv = TRUE)

# Ajustar el modelo PGLS
pgls_model <- pgls(log_svl ~ real_reg, data = comparative_data)
summary(pgls_model)

# Obtener los residuos del modelo PGLS
res_svl <- residuals(pgls_model)
dt74$adjusted_body_size <- residuals(pgls_model)

# Realizar un boxplot de los valores ajustados
ggplot(data, aes(x = clade, y = adjusted_body_size)) +
  geom_boxplot() +
  theme_minimal()
