###################################
#   Preparing BIOMES
#################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/My analyses/1. Distributions")
library(raster)
library(ncdf4)
library(rgdal)

# 1. Import distributions
#=========================
distAM <- readOGR('dist_America/dist_America.shp')
distEA <- readOGR('dist_Eurasia/dist_Eurasia.shp')

# 2. Import Ecoregions
#=====================
library(rgeos)
ecoregs <- readOGR('Ecoregions2017/Ecoregions2017.shp')


# 3. Extract biomes by distributions
#========================================

 # EURASIA RADIATION
 #==================
overlap1 <- intersect(distEA, ecoregs)
plot(overlap1)
str(overlap1, max.level=2)

library(gdata)
overlap1@data <- drop.levels(overlap1@data)

str(overlap1@data)
str(overlap1@polygons, max.level=2)

# Extract areas
areas <- sapply(slot(overlap1, "polygons"), slot, "area")

# Create an empty dataframe of biomes
biome.per.lineage <- data.frame(biome=levels(as.factor(overlap1@data$BIOME_NAME)))

a = 1
for(i in levels(as.factor(overlap1@data$sci_name))){
  
  biomes.in.poly <- overlap1@data$BIOME_NAME[overlap1@data$sci_name == i]
  areas.poly <- areas[overlap1@data$sci_name == i]
  
  value <- tapply(areas.poly, biomes.in.poly, sum)
  
  percentage <- value / sum(value, na.rm=T) * 100
  
  index <- match(names(percentage), biome.per.lineage$biome)
  biome.per.lineage[index, "pcol"] <- percentage
  
  colnames(biome.per.lineage)[a+1] <- i
  
  a = a + 1
}

# Customize dataframe
rownames(biome.per.lineage) <- biome.per.lineage$biome
biome.per.lineage <- biome.per.lineage[,-1]
biome.per.lineage <- data.frame(t(biome.per.lineage))
biome.per.lineage[is.na(biome.per.lineage)] <- 0
biome.per.lineage <- round(biome.per.lineage, 3)
biome.per.lineage$sp_lin <- rownames(biome.per.lineage)
biome.per.lineage <- biome.per.lineage[,c(ncol(biome.per.lineage),1:(ncol(biome.per.lineage)-1))] # to be generizable
rownames(biome.per.lineage) <- NULL


# Create a new colum of the final biome for each species 
biome.per.lineage$s.biom <- NA
for(i in 1:nrow(biome.per.lineage)){
  bms <- ncol(biome.per.lineage) - 1
  biome.per.lineage$s.biom[i] <- colnames(biome.per.lineage[i,2:bms])[which(biome.per.lineage[i,2:bms] == max(biome.per.lineage[i,2:bms]))]
}

biome.per.lineage$sp_lin
table(biome.per.lineage$s.biom)

# SAVE
# write.csv(biome.per.lineage, file='biome_EA.csv', row.names=F)


# REPEAT FOR AMERICAN RADIATION
#================================
overlap2 <- intersect(distAM, ecoregs)
plot(overlap2)
str(overlap2, max.level=2)

library(gdata)
overlap2@data <- drop.levels(overlap2@data)

str(overlap2@data)
str(overlap2@polygons, max.level=2)

# Extract areas
areas <- sapply(slot(overlap2, "polygons"), slot, "area")

# Create an empty dataframe of biomes
biome.per.lineage <- data.frame(biome=levels(as.factor(overlap2@data$BIOME_NAME)))

a = 1
for(i in levels(as.factor(overlap2@data$sci_name))){
  
  biomes.in.poly <- overlap2@data$BIOME_NAME[overlap2@data$sci_name == i]
  areas.poly <- areas[overlap2@data$sci_name == i]
  
  value <- tapply(areas.poly, biomes.in.poly, sum)
  
  percentage <- value / sum(value, na.rm=T) * 100
  
  index <- match(names(percentage), biome.per.lineage$biome)
  biome.per.lineage[index, "pcol"] <- percentage
  
  colnames(biome.per.lineage)[a+1] <- i
  
  a = a + 1
}

# Customize dataframe
rownames(biome.per.lineage) <- biome.per.lineage$biome
biome.per.lineage <- biome.per.lineage[,-1]
biome.per.lineage <- data.frame(t(biome.per.lineage))
biome.per.lineage[is.na(biome.per.lineage)] <- 0
biome.per.lineage <- round(biome.per.lineage, 3)
biome.per.lineage$sp_lin <- rownames(biome.per.lineage)
biome.per.lineage <- biome.per.lineage[,c(ncol(biome.per.lineage),1:(ncol(biome.per.lineage)-1))] # to be generizable
rownames(biome.per.lineage) <- NULL


# Create a new colum of the final biome for each species 
biome.per.lineage$s.biom <- NA
for(i in 1:nrow(biome.per.lineage)){
  bms <- ncol(biome.per.lineage) - 1
  biome.per.lineage$s.biom[i] <- colnames(biome.per.lineage[i,2:bms])[which(biome.per.lineage[i,2:bms] == max(biome.per.lineage[i,2:bms]))]
}

biome.per.lineage$sp_lin
table(biome.per.lineage$s.biom)

# SAVE
# write.csv(biome.per.lineage, file='biome_AM.csv', row.names=F)



# 4. UNIQUE BIOMES (Group into a few categories)
#================================================
biomeEA <- read.csv("biome_EA.csv", sep=";")
biomeAM <- read.csv("biome_AM.csv", sep=";")

biomeEA <- biomeEA[,c(1,14:15)]
biomeAM <- biomeAM[,c(1,17:18)]

# Merge datasets
colnames(biomeEA) <- colnames(biomeAM)
biome <- rbind(biomeEA, biomeAM)
unique(biome$general_biome)

# Regroup in 5 categories
# 1 = Boreal Forest
# 2 = Desert
# 3 = Grasslands
# 4 = Temperate forest
# 5 = Tropical forest

biome$new_biome <- NA
for (i in 1:nrow(biome)){
  if (biome[i,3]=="Mediterranean forest" |
      biome[i,3]=="Temperate conifer forest" |
      biome[i,3]=="Temperate broadleaf forest"){
    biome$new_biome[i] <- "Temperate Forest"
  } else if (biome[i,3]=="Tropical moist forest" |
             biome[i,3]=="Tropical coniferous forest" |
             biome[i,3]=="Tropical dry forest"){
    biome$new_biome[i] <- "Tropical forest"
  } else if (biome[i,3]=="Montane grasslands" | biome[i,3]=="Temperate grassland"){
    biome$new_biome[i] <- "Grasslands"
  } else {
    biome$new_biome[i] <- biome$general_biome[i]
  }
}

unique(biome$new_biome)
biome$sp_lin  # 75 sp

# SAVE
# write.csv(biome, "biome_all.csv")


# 5. Import phylogeny
#=====================
library(phytools)
library(ape)
tr <- read.tree("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/My analyses/PhyloRana_75sp.tre")

# Make sp match
rownames(biome) <- biome$sp_lin
tr$tip.label == biome$sp_lin
biome <- biome[tr$tip.label,]

b <- setNames(biome$new_biome, rownames(biome))
dotTree(tr, b)

bmode <- as.factor(b)
levels(bmode)
contMap(tr, bmode, type = "fan", fsize=0.7)

