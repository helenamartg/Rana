#############################
#   Cleaning morpho Ranidae
#############################
rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/My analyses")

# 1. Import data
#=================
svl <- read.csv("3. Traits_dataset/traitsRana_HMG.csv", sep=";")
svl <- svl[c(1:77), c(1:5,7,9)]
mor <- read.csv("3. Traits_dataset/Morpho_rana.csv", sep=";")
mor <- mor[c(1:77),]

# 2. Combine datasets
#=====================
dt <- NULL
for (i in 1:nrow(mor)){
  for (j in 1:nrow(svl)){
    if (mor[i,1]==svl[j,1]){
      x <- cbind(svl[j,c(1:5)], mor[i,], svl[j,c(6:7)])
      dt <- rbind(dt, x)
      
    }
  }
}
dt

# 3. Obtain a dataset  for all species (with SVL, clutch size, longevity and ecotype)
#======================================
# Remove repeated columns
dt <- dt[,-c(6,7)]
unique(dt$ECOTYPE)
dt_red <- dt[,-c(6:26)]

# Remove Rana kunyensis (= R. coreana)
dt_red <- dt_red[dt_red$Species!="Rana_kunyuensis",]

# Remove Rana tlaloci (no climatic data)
dt_red <- dt_red[dt_red$Species!="Rana_tlaloci",]

# Remove Rana psilonota (no SVL data)
dt_red <- dt_red[dt_red$Species!="Rana_psilonota",]

dt74 <- dt_red

# Add the real region
colnames(dt74)
colnames(dt74) <- c("sp", "genus", "phylo_reg", "state", "SVL", "clutch_size", "ecotype")

dt74$real_reg <- NA
for (i in 1:nrow(dt74)){
  if (dt74[i,4]=="AS" | dt74[i,4]=="EU"){
    dt74[i,8] <- "EA"
  } else {
    dt74[i,8] <- "AM"
  }
}

# Reorder
dt74 <- dt74[, c(1:4,8,7,6,5)]
dim(dt74)

# SAVE the reduced dataset with 74sp 
# write.csv(dt74, "3. Traits_dataset/dt_74sp.csv")


# 4. Obtain a dataset of morphological variables
#=================================================
morpho <- dt
rownames(morpho) <- dt$Species

# Missing values
apply(morpho, 2, function(x){sum(is.na(x))})

# transpose matrix to have sp names as columns and count the missing values per species
transp <- as.matrix(t(morpho))
colnames(transp) <- dt$Species
transp <- as.data.frame(transp)
apply(transp, 2, function(x){sum(is.na(x))})

# remove species that have more than 20 missing values because are sp with any measure
new.morpho <- NULL
for (i in 1:ncol(transp)){
    if (sum(is.na(transp[,i])) < 20){
      x <- transp[,i]
      new.morpho <- cbind(new.morpho, x)
  }
}

colnames(new.morpho) <- new.morpho[1,]
new.morpho <- new.morpho[-1,]
new.morpho <- as.data.frame(new.morpho)

# Transpose again 
morpho2 <- as.matrix(t(new.morpho))
colnames(morpho2) <- colnames(morpho[-1]) 
dim(morpho2)
# The new dataset morpho2 has 55 species.

# See how many NAs are in each variable
apply(morpho2, 2, function(x){sum(is.na(x))}) 

# Leave only morphological variables with less than 10 NAs: HL, HW, FL, TL
def_morph <- morpho2[,c(1:4, 6:7, 20, 22, 26, 27)]
def_morph <- as.data.frame(def_morph)

apply(def_morph, 2, function(x){sum(is.na(x))}) # HW for all the 55 sp

# Reorder
def_morph$sp <- rownames(def_morph)
dt55 <- def_morph[,c(11,1:3,9,10,4:8)]
colnames(dt55) <- c("sp", "genus", "phylo_reg", "state", "Clutch_size", "ecotype",
                    "SVL", "HL", "HW", "FL", "TL")

# SAVE this dataset (ONly HW complete)
# write.csv(dt55, '3. Traits_dataset/dt_morpho_55sp.csv')


# 5. Exploratory plots
#======================
dt55 <- dt55[,-1]

# Standarize and log all variables
dt74$log.SVL <- as.vector(log(dt74$SVL))
dt55[,8:10] <- log(dt55[,8:10])
dt55$SVL <- log(dt55$SVL)

# Residuals only for HW
rownames(dt55) <- dt55$sp
dt55$res_HW <- resid(lm(dt55$HW ~ dt55$SVL))

# Plot SVL by region in one single plot
par(mar=c(4,4,4,4))
boxplot(dt74$SVL~ dt74$real_reg, col="white")
stripchart(dt74$SVL~ dt74$real_reg, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = c("gold", "darkorchid"))

# Plot SVL by state
pal <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
boxplot(dt74$SVL ~ dt74$state, col=pal)

# Plot togheter
par(mfrow=c(1,2), mar=c(4,4,4,2))
boxplot(dt74$SVL~ dt74$real_reg, col="white", xlab="", ylab="SVL (mm)")
stripchart(dt74$SVL~ dt74$real_reg, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = c("gold", "darkorchid"))
boxplot(dt74$SVL ~ dt74$state, col=pal, xlab="", ylab="")


# Plot the SVL by ecotype
layout(1,1)
unique(dt74$ecotype)
pal2 <- c("dodgerblue4", "cornflowerblue", "chocolate4")
boxplot(dt74$SVL~ dt74$ecotype, col="white", xlab="", ylab="SVL (mm)", lab.cex=0.5,
        main="76 sp")
stripchart(dt74$SVL~ dt74$ecotype, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col=pal2)

# Plot HW by ecotype 
boxplot(dt55$res_HW ~ dt55$ecotype, col="white", xlab="", ylab="HW (mm)", lab.cex=0.5, main="55 sp")
stripchart(dt55$res_HW~ dt55$ecotype, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col=pal2)


# Import phylogeny
library(phytools)
library(ape)
phylorana <- read.tree("PhyloRana_74sp.tre")
phylorana$tip.label
phylorana <- drop.tip(phylorana, 55)
rownames(dt74) <- dt74$sp
svl <- setNames(dt74$SVL, rownames(dt74))

# Make dataset match with phylogeny
rownames(dt74) ==  phylorana$tip.label
dt74 <- dt74[phylorana$tip.label,]

#plot svl across tree
contMap(phylorana, log(svl), fsize=.5, type='fan')

layout(1,1)
obj<-contMap(phylorana, svl, plot=FALSE)
plotTree.wBars(obj$tree, svl, fsize=0.6, scale=0.2,tip.labels=TRUE,
               method="plotSimmap",colors=obj$cols)



