#############################
#   Cleaning morpho Ranidae
#############################
rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

# 1. Import data
#=================
svl <- read.csv("DATA/Raw_data/traitsRana_HMG.csv", sep=";")
svl <- svl[c(1:77), c(1:5,7,9)]
mor <- read.csv("DATA/Raw_data/Morpho_rana.csv", sep=";")
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

# Log-transformed body size and clutch
dt74$log_svl <- log(dt74$SVL)
dt74$log_clutch <- log(dt74$clutch_size)

# Residuals of clutch size
dt74$resid_clutch <- resid(lm(dt74$log_clutch ~ dt74$log_svl, na.action = na.exclude))

# SAVE the reduced dataset with 74sp 
# write.csv(dt74, "DATA/dt_74sp.csv")


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

# log-transformed
dt55 <- as.data.frame(dt55)
dt55$log_svl <- log(as.numeric(dt55$SVL))

# calculate residuals for all morphological variables
only.morpho <- dt55[,8:11]
res.morpho <- matrix(nrow=nrow(dt55), ncol=ncol(only.morpho))
for (i in 1:ncol(only.morpho)){
  res.morpho[,i] <- resid(lm(log(as.numeric(only.morpho[,i])) ~ dt55$log_svl, na.action = na.exclude))
}
colnames(res.morpho) <- paste("res", colnames(only.morpho), sep = "_")
rownames(res.morpho) <- rownames(only.morpho)

# Merge datasets
merged.morph <- cbind(dt55, res.morpho)

# Residuals of clutch size
merged.morph$res_clutch <- resid(lm(log(as.numeric(merged.morph$Clutch_size)) ~ merged.morph$log_svl, na.action = na.exclude))

# Reorder
final.morph <- merged.morph[,c(1:4,6,5,7:17)]

# SAVE this dataset (ONly HW complete)
# write.csv(final.morph, 'DATA/dt_morpho_55sp.csv')


