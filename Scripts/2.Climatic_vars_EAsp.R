####################################
#     Getting climatic variables
#         Eurasian radiation
####################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/My analyses/")

library(maptools)
library(raster)
library(ncdf4)
library(rgdal)

# 1. Import distributions
#========================
distEA <- readOGR('1. Distributions/dist_Eurasia/dist_Eurasia.shp')
distAM <- readOGR('1. Distributions/dist_America/dist_America.shp')

# proj4string(distEA) <- CRS("+init=epsg:4326")


# 2. Load climatic rasters
#==========================
# NDVI
NDVI_raster <-'2. Climatic_variables/BIOCLIM/MOD_NDVI_M_2023-01-01_rgb_1440x720.TIFF' 
ndvi <- raster(NDVI_raster)
library(RColorBrewer)
plot(ndvi, main='primariy productivity (NDVI)', col=brewer.pal(9, name="YlGn"))

# Aridity index
aridity_raster <- '2. Climatic_variables/BIOCLIM/Global-AI_ET0_v3_annual/AI_RES025_correct.tif'
arid <- raster(aridity_raster)
plot(arid, main='aridity index', col=brewer.pal(9, name="RdYlBu"))
plot(distEA, add=T)


# Potential evapotranspiration
PET_raster <- '2. Climatic_variables/BIOCLIM/Global-AI_ET0_v3_annual/ET0_RES025.tif'
pet <- raster(PET_raster)
plot(pet, main='potential evapotranspiration', col=topo.colors(10))


# Elevation (WorldClim)
elevation <- '2. Climatic_variables/BIOCLIM/wc2.1_10m_elev.tif'
elev <- raster(elevation)
plot(elevation, main='Elevation')
plot(distAM, add=T)
plot(distEA, add=T)


# prec, srad, tavg, tmax, tmin, vapr, wind
prec.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_prec', 
                     pattern='.tif', full.names=T)
prec <- stack(prec.f)

tmax.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_tmax', pattern='.tif', full.names=T)
tmax <- stack(tmax.f)

vapr.f <- list.files('2. Climatic_variables/BIOCLIM//wc2.1_10m_vapr', pattern='.tif', full.names=T)
vapr <- stack(vapr.f)

tmin.f <- list.files('2. Climatic_variables/BIOCLIM//wc2.1_10m_tmin', pattern='.tif', full.names=T)
tmin <- stack(tmin.f)

tavg.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_tavg', pattern='.tif', full.names=T)
tavg <- stack(tavg.f)

wind.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_wind', pattern='.tif', full.names=T)
wind <- stack(wind.f)

srad.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_srad', pattern='.tif', full.names=T)
srad <- stack(srad.f)

plot(tmax[[1]])


# 3. EXTRACT VALUES PER SP
#==========================


#######
# NDVI
#######
ndvi.vals <- extract(ndvi, distEA, na.rm=T)
names(ndvi.vals) <- distEA@data$sci_name
str(ndvi.vals)

ndvi.save <- ndvi.vals

NDVIavg <- data.frame(matrix(ncol=1, nrow=length(ndvi.vals)))
colnames(NDVIavg) <- 'NDVIavg'
NDVImax <- data.frame(matrix(ncol=1, nrow=length(ndvi.vals)))
colnames(NDVImax) <- 'NDVImax'
NDVImin <- data.frame(matrix(ncol=1, nrow=length(ndvi.vals)))
colnames(NDVImin) <- 'NDVImin'

NDVIsd <- data.frame(matrix(ncol=1, nrow=length(ndvi.vals)))
colnames(NDVIsd) <- 'NDVIsd'

for(i in 1:length(ndvi.vals)){
  NDVIavg[i,] <- mean(ndvi.vals[[i]], na.rm=T)
  NDVImax[i,] <- max(ndvi.vals[[i]], na.rm=T)
  NDVImin[i,] <- min(ndvi.vals[[i]], na.rm=T)
  
  NDVIsd[i,] <- sd(ndvi.vals[[i]], na.rm=T)
}


##########
# Aridity
##########
arid.vals <- extract(arid, distEA, na.rm=T)
names(arid.vals) <- distEA@data$sci_name
str(arid.vals)

arid.save <- arid.vals

ARIDavg <- data.frame(matrix(ncol=1, nrow=length(arid.vals)))
colnames(ARIDavg) <- 'ARIDavg'
ARIDmax <- data.frame(matrix(ncol=1, nrow=length(arid.vals)))
colnames(ARIDmax) <- 'ARIDmax'
ARIDmin <- data.frame(matrix(ncol=1, nrow=length(arid.vals)))
colnames(ARIDmin) <- 'ARIDmin'

aridSD <- data.frame(matrix(ncol=1, nrow=length(arid.vals)))
colnames(aridSD) <- 'ARIDsd'

for(i in 1:length(arid.vals)){
  ARIDavg[i,] <- mean(arid.vals[[i]], na.rm=T)
  ARIDmax[i,] <- max(arid.vals[[i]], na.rm=T)
  ARIDmin[i,] <- min(arid.vals[[i]], na.rm=T)
  
  aridSD[i,] <- sd(arid.vals[[i]], na.rm=T)
}


##########
# PET
##########
pet.vals <- extract(pet, distEA, na.rm=T)
names(pet.vals) <- distEA@data$sci_name
str(pet.vals)

pet.save <- pet.vals

PETavg <- data.frame(matrix(ncol=1, nrow=length(pet.vals)))
colnames(PETavg) <- 'PETavg'
PETmax <- data.frame(matrix(ncol=1, nrow=length(pet.vals)))
colnames(PETmax) <- 'PETmax'
PETmin <- data.frame(matrix(ncol=1, nrow=length(pet.vals)))
colnames(PETmin) <- 'PETmin'

petSD <- data.frame(matrix(ncol=1, nrow=length(pet.vals)))
colnames(petSD) <- 'PETsd'

for(i in 1:length(pet.vals)){
  PETavg[i,] <- mean(pet.vals[[i]], na.rm=T)
  PETmax[i,] <- max(pet.vals[[i]], na.rm=T)
  PETmin[i,] <- min(pet.vals[[i]], na.rm=T)
  
  petSD[i,] <- sd(pet.vals[[i]], na.rm=T)
}


##########
# elev
#########
elev.vals <- extract(elev, distEA, na.rm=T)
names(elev.vals) <- distEA@data$sci_name
str(elev.vals)

elev.save <- elev.vals
elev.vals <- elev.save

ELEVavg <- data.frame(matrix(ncol=1, nrow=length(elev.vals)))
colnames(ELEVavg) <- 'ELEVavg'
ELEVmax <- data.frame(matrix(ncol=1, nrow=length(elev.vals)))
colnames(ELEVmax) <- 'ELEVmax'
ELEVmin <- data.frame(matrix(ncol=1, nrow=length(elev.vals)))
colnames(ELEVmin) <- 'ELEVmin'

ELEVsd <- data.frame(matrix(ncol=1, nrow=length(elev.vals)))
colnames(ELEVsd) <- 'ELEVsd'

for(i in 1:length(elev.vals)){
  ELEVavg[i,] <- mean(elev.vals[[i]], na.rm=T)
  ELEVmax[i,] <- max(elev.vals[[i]], na.rm=T)
  ELEVmin[i,] <- min(elev.vals[[i]], na.rm=T)
  
  ELEVsd[i,] <- sd(elev.vals[[i]], na.rm=T)
}

#########
# tmax
#########
tmax.vals <- extract(tmax, distEA, na.rm=T)
names(tmax.vals) <- distEA@data$sci_name
str(tmax.vals)

tmax.save <- tmax.vals

tmax.mean <- data.frame(matrix(ncol=12, nrow=length(tmax.vals)))
tmax.sd <- data.frame(matrix(ncol=12, nrow=length(tmax.vals)))
for(i in 1:length(tmax.vals)){
  tmax.mean[i,] <- apply(tmax.vals[[i]], 2, mean, na.rm=T)
  tmax.sd[i,] <- apply(tmax.vals[[i]], 2, sd, na.rm=T)
}
tmax.mean <- t(tmax.mean)
tmax.sd <- t(tmax.sd)

tmax.mean.annual <- apply(tmax.mean, 2, mean)
tmax.sd.annual <- apply(tmax.mean, 2, sd)

TMAXavg <- tmax.mean.annual
TMAXsd <- tmax.sd.annual


#########
# tmin
#########
tmin.vals <- extract(tmin, distEA, na.rm=T)
names(tmin.vals) <- distEA@data$sci_name
str(tmin.vals)

tmin.save <- tmin.vals

tmin.mean <- data.frame(matrix(ncol=12, nrow=length(tmin.vals)))
tmin.sd <- data.frame(matrix(ncol=12, nrow=length(tmin.vals)))
for(i in 1:length(tmin.vals)){
  tmin.mean[i,] <- apply(tmin.vals[[i]], 2, mean, na.rm=T)
  tmin.sd[i,] <- apply(tmin.vals[[i]], 2, sd, na.rm=T)
}
tmin.mean <- t(tmin.mean)
tmin.sd <- t(tmin.sd)

tmin.mean.annual <- apply(tmin.mean, 2, mean)
tmin.sd.annual <- apply(tmin.mean, 2, sd)

TMINavg <- tmin.mean.annual
TMINsd <- tmin.sd.annual


########
# TAVG
########
tavg.vals <- extract(tavg, distEA, na.rm=T)
names(tavg.vals) <- distEA@data$sci_name
str(tavg.vals)

tavg.save <- tavg.vals

tavg.mean <- data.frame(matrix(ncol=12, nrow=length(tavg.vals)))
tavg.sd <- data.frame(matrix(ncol=12, nrow=length(tavg.vals)))
for(i in 1:length(tavg.vals)){
  tavg.mean[i,] <- apply(tavg.vals[[i]], 2, mean, na.rm=T)
  tavg.sd[i,] <- apply(tavg.vals[[i]], 2, sd, na.rm=T)
}
tavg.mean <- t(tavg.mean)
tavg.sd <- t(tavg.sd)

tavg.mean.annual <- apply(tavg.mean, 2, mean)
tavg.sd.annual <- apply(tavg.mean, 2, sd)

TAVGavg <- tavg.mean.annual
TAVGsd <- tavg.sd.annual


########
# VAPR
########
vapr.vals <- extract(vapr, distEA, na.rm=T)
names(vapr.vals) <- distEA@data$sci_name
str(vapr.vals)

vapr.save <- vapr.vals

vapr.mean <- data.frame(matrix(ncol=12, nrow=length(vapr.vals)))
vapr.sd <- data.frame(matrix(ncol=12, nrow=length(vapr.vals)))
for(i in 1:length(vapr.vals)){
  vapr.mean[i,] <- apply(vapr.vals[[i]], 2, mean, na.rm=T)
  vapr.sd[i,] <- apply(vapr.vals[[i]], 2, sd, na.rm=T)
}
vapr.mean <- t(vapr.mean)
vapr.sd <- t(vapr.sd)

vapr.mean.annual <- apply(vapr.mean, 2, mean)
vapr.sd.annual <- apply(vapr.mean, 2, sd)
VAPRmax <- apply(vapr.mean, 2, max)
VAPRmin <- apply(vapr.mean, 2, min)

VAPRavg <- vapr.mean.annual
VAPRsd <- vapr.sd.annual


###########
# WIND
##########
wind.vals <- extract(wind, distEA, na.rm=T)
names(wind.vals) <- distEA@data$sci_name
str(wind.vals)

wind.save <- wind.vals

wind.mean <- data.frame(matrix(ncol=12, nrow=length(wind.vals)))
wind.sd <- data.frame(matrix(ncol=12, nrow=length(wind.vals)))
for(i in 1:length(wind.vals)){
  wind.mean[i,] <- apply(wind.vals[[i]], 2, mean, na.rm=T)
  wind.sd[i,] <- apply(wind.vals[[i]], 2, sd, na.rm=T)
}
wind.mean <- t(wind.mean)
wind.sd <- t(wind.sd)

wind.mean.annual <- apply(wind.mean, 2, mean)
wind.sd.annual <- apply(wind.mean, 2, sd)
WINDmax <- apply(wind.mean, 2, max)
WINDmin <- apply(wind.mean, 2, min)

WINDavg <- wind.mean.annual
WINDsd <- wind.sd.annual


# IMPORTANT
# NOTE: for precipitation and solar radiation we will take a different approach:
# first we will add all the layers to estimate the anual precipitation and s. radiation

###########
# PREC
##########
# first sum all the layers
sum.prec <- sum(prec)

prec.vals <- extract(sum.prec, distEA, na.rm=T)
names(prec.vals) <- distEA@data$sci_name
str(prec.vals)

prec.save <- prec.vals

PRECavg <- data.frame(matrix(ncol=1, nrow=length(prec.vals)))
colnames(PRECavg) <- 'PRECavg'
PRECmax <- data.frame(matrix(ncol=1, nrow=length(prec.vals)))
colnames(PRECmax) <- 'PRECmax'
PRECmin <- data.frame(matrix(ncol=1, nrow=length(prec.vals)))
colnames(PRECmin) <- 'PRECmin'

PRECsd <- data.frame(matrix(ncol=1, nrow=length(prec.vals)))
colnames(PRECsd) <- 'PRECsd'

for(i in 1:length(prec.vals)){
  PRECavg[i,] <- mean(prec.vals[[i]], na.rm=T)
  PRECmax[i,] <- max(prec.vals[[i]], na.rm=T)
  PRECmin[i,] <- min(prec.vals[[i]], na.rm=T)
  
  PRECsd[i,] <- sd(prec.vals[[i]], na.rm=T)
}


#####################
# srad
#####################

sum.srad <- sum(srad)

srad.vals <- extract(sum.srad, distEA, na.rm=T)
names(srad.vals) <- distEA@data$sci_name
str(srad.vals)

srad.save <- srad.vals

SRADavg <- data.frame(matrix(ncol=1, nrow=length(srad.vals)))
colnames(SRADavg) <- 'SRADavg'
SRADmax <- data.frame(matrix(ncol=1, nrow=length(srad.vals)))
colnames(SRADmax) <- 'SRADmax'
SRADmin <- data.frame(matrix(ncol=1, nrow=length(srad.vals)))
colnames(SRADmin) <- 'SRADmin'

SRADsd <- data.frame(matrix(ncol=1, nrow=length(srad.vals)))
colnames(SRADsd) <- 'SRADsd'

for(i in 1:length(srad.vals)){
  SRADavg[i,] <- mean(srad.vals[[i]], na.rm=T)
  SRADmax[i,] <- max(srad.vals[[i]], na.rm=T)
  SRADmin[i,] <- min(srad.vals[[i]], na.rm=T)
  
  SRADsd[i,] <- sd(srad.vals[[i]], na.rm=T)
}



##########################################
# Organize and save the data
##########################################

env95 <- as.data.frame(cbind(ARIDavg, ARIDmax, ARIDmin, aridSD,
                             NDVIavg, NDVImax, NDVImin, NDVIsd,
                             PETavg, PETmax, PETmin, petSD,
                             PRECavg, PRECmax, PRECmin, PRECsd, 
                             SRADavg, SRADmax, SRADmin, SRADsd,
                             TAVGavg, TAVGsd, 
                             TMAXavg, TMAXsd, 
                             TMINavg, TMINsd, 
                             VAPRavg, VAPRmax, VAPRmin, VAPRsd,
                             ELEVavg, ELEVmax, ELEVmin, ELEVsd,
                             WINDavg, WINDmax, WINDmin, WINDsd,
                             SRADavg, SRADmax, SRADmin, SRADsd))
colnames(env95)
rownames(env95) <- distEA@data$sci_name
rownames(env95)

# write.csv(env95, "env_vars_EA_sp.csv")


