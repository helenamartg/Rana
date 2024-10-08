rm(list = ls())

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

###################
# Preparing NDVI 
###################
library(raster)
library(terra)
library(dplyr)

# load distributions (one polygon per species)
dists_am <- vect("DATA/dist_America/dist_America.shp")

#NDVI annual
NDVI.f <- list.files("DATA/NDVI_GeoTIFFs",
                     pattern=".TIFF", full.name=T) #all months for 23 years  (276)
NDVI <- rast(NDVI.f) #stack all of them
class(NDVI)
class(dists_am)

# extract all values within ranges by lineage
ndvi.vals <- extract(NDVI, dists_am, na.rm=T)

# make a copy of extracted values (just in case and to make sure that everything goes ok)
ndvi.save <- ndvi.vals
ndvi.vals <- ndvi.save

#mean by groups - ID 
result <- aggregate(. ~ ID, data = ndvi.save, FUN = function(x) mean(x, na.rm = TRUE))
ndvi.mean <- result

unique(result$ID)

sp.names <- values(dists_am)
sp.names <- sp.names$sci_name

rownames(ndvi.mean) <- sp.names
ndvi.mean <- ndvi.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
ndvi.avg <- apply(ndvi.mean, 1, mean, na.rm=T)
ndviMmax <- apply(ndvi.mean, 1, max, na.rm=T)
ndviMmin <- apply(ndvi.mean, 1, min, na.rm=T)
ndviSD <- apply(ndvi.mean, 1, sd, na.rm=T)

ndvi_total <- data.frame(sp.names,ndvi.avg, ndviMmax, ndviMmin, ndviSD)
# write.csv(ndvi_total, "DATA/ndvi_annual_AM.csv")



# Same for Eurasian species
#===========================
# load distributions (one polygon per species)
dists_ea <- vect("DATA/dist_Eurasia/dist_Eurasia.shp")

#NDVI annual
NDVI.f <- list.files("DATA/NDVI_GeoTIFFs",
                     pattern=".TIFF", full.name=T) #all months for 23 years  (276)
NDVI <- rast(NDVI.f) #stack all of them
class(NDVI)
class(dists_ea)

# extract all values within ranges by lineage
ndvi.vals <- extract(NDVI, dists_ea, na.rm=T)

# make a copy of extracted values (just in case and to make sure that everything goes ok)
ndvi.save <- ndvi.vals
ndvi.vals <- ndvi.save

#mean by groups - ID 
result <- aggregate(. ~ ID, data = ndvi.save, FUN = function(x) mean(x, na.rm = TRUE))
ndvi.mean <- result

unique(result$ID)

sp.names <- values(dists_ea)
sp.names <- sp.names$sci_name

rownames(ndvi.mean) <- sp.names
ndvi.mean <- ndvi.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
ndvi.avg <- apply(ndvi.mean, 1, mean, na.rm=T)
ndviMmax <- apply(ndvi.mean, 1, max, na.rm=T)
ndviMmin <- apply(ndvi.mean, 1, min, na.rm=T)
ndviSD <- apply(ndvi.mean, 1, sd, na.rm=T)

ndvi_total <- data.frame(sp.names,ndvi.avg, ndviMmax, ndviMmin, ndviSD)
# write.csv(ndvi_total, "DATA/ndvi_annual_EA.csv")


##############################################################################
##### MERGE with climatic dataset
###############################################################################
clim.am <- read.csv("DATA/clim_AM_37sp.csv")
rownames(clim.am) <- clim.am$X

clim.ea <- read.csv("DATA/clim_EA_31sp.csv")
rownames(clim.ea) <- clim.ea$X

ndvi.am <- read.csv("DATA/ndvi_annual_AM.csv")
rownames(ndvi.am) <- ndvi.am$sp.names

ndvi.ea <- read.csv("DATA/ndvi_annual_EA.csv")
rownames(ndvi.ea) <- ndvi.ea$sp.names

# Check if sp names match
rownames(clim.am)==rownames(ndvi.am)
rownames(clim.ea)==rownames(ndvi.ea)
ndvi.ea <- ndvi.ea[clim.ea$X,]

# sustituir antiguos valores de ndvi por nuevos
clim.am[,7:10] <- ndvi.am[,3:6]
clim.ea[,7:10] <- ndvi.ea[,3:6]



#SAVE
# write.csv(clim.am, "DATA/clim_AM_37sp.csv")
# write.csv(clim.ea, "DATA/clim_EA_31sp.csv")


