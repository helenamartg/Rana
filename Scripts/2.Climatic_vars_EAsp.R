####################################
#     Getting climatic variables
#         Eurasian radiation
####################################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/My analyses/")

library(raster)
library(ncdf4)
library(terra)
library(RColorBrewer)

# 1. Import distributions
#========================
distEA <- vect('1. Distributions/dist_Eurasia/dist_Eurasia.shp')

# 2. Load climatic rasters
#==========================
# NDVI
NDVI <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/ndvi_annual_EA.csv")
rownames(NDVI) <- NDVI$sp.names

# Aridity index
aridity_raster <- '2. Climatic_variables/BIOCLIM/Global-AI_ET0_v3_annual/AI_RES01.tif'
arid <- rast(aridity_raster)
plot(arid, main='aridity index', col=brewer.pal(9, name="RdYlBu"))
plot(distEA, add=T)

# Potential evapotranspiration
PET_raster <- '2. Climatic_variables/BIOCLIM/Global-AI_ET0_v3_annual/ET0_RES01.tif'
pet <- rast(PET_raster)
plot(pet, main='potential evapotranspiration', col=topo.colors(10))


# Elevation (WorldClim)
elevation <- '2. Climatic_variables/BIOCLIM/wc2.1_10m_elev.tif'
elev <- rast(elevation)
plot(elevation, main='Elevation')
plot(distEA, add=T)
plot(distEA, add=T)


# prec, srad, tavg, tmax, tmin, vapr, wind
prec.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_prec', 
                     pattern='.tif', full.names=T)
prec <- rast(prec.f)

tmax.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_tmax', pattern='.tif', full.names=T)
tmax <- rast(tmax.f)

vapr.f <- list.files('2. Climatic_variables/BIOCLIM//wc2.1_10m_vapr', pattern='.tif', full.names=T)
vapr <- rast(vapr.f)

tmin.f <- list.files('2. Climatic_variables/BIOCLIM//wc2.1_10m_tmin', pattern='.tif', full.names=T)
tmin <- rast(tmin.f)

tavg.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_tavg', pattern='.tif', full.names=T)
tavg <- rast(tavg.f)

wind.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_wind', pattern='.tif', full.names=T)
wind <- rast(wind.f)

srad.f <- list.files('2. Climatic_variables/BIOCLIM/wc2.1_10m_srad', pattern='.tif', full.names=T)
srad <- rast(srad.f)

plot(tmax[[1]])


# 3. EXTRACT VALUES PER SP
#==========================
#######
# NDVI
#######
# Done in previous Script when downloading NDVI rasters from NASA webpage
NDVIavg <- NDVI$ndvi.avg
NDVImax <- NDVI$ndviMmax
NDVImin <- NDVI$ndviMmin
ndviSD <- NDVI$ndviSD

##########
# Aridity
##########
arid.vals <- extract(arid, distEA, na.rm=T)
str(arid.vals)

arid.save <- arid.vals

#mean by groups - ID 
ARIDavg <- aggregate(. ~ ID, data = arid.save, FUN = function(x) mean(x, na.rm = TRUE))
ARIDmin <- aggregate(. ~ ID, data = arid.save, FUN = function(x) min(x, na.rm = TRUE))
ARIDmax <- aggregate(. ~ ID, data = arid.save, FUN = function(x) max(x, na.rm = TRUE))
aridSD <- aggregate(. ~ ID, data = arid.save, FUN = function(x) sd(x, na.rm = TRUE))

unique(ARIDavg$ID)
unique(aridSD$ID)

sp.names <- values(distEA)
sp.names <- sp.names$sci_name

ARID <- cbind(ARIDavg = ARIDavg[,2], ARIDmin = ARIDmin[,2], ARIDmax = ARIDmax[,2], aridSD = aridSD[,2])
rownames(ARID) <- sp.names


##########
# PET
##########
pet.vals <- extract(pet, distEA, na.rm=T)
str(pet.vals)

pet.save <- pet.vals

#mean by groups - ID 
PETavg <- aggregate(. ~ ID, data = pet.save, FUN = function(x) mean(x, na.rm = TRUE))
PETmin <- aggregate(. ~ ID, data = pet.save, FUN = function(x) min(x, na.rm = TRUE))
PETmax <- aggregate(. ~ ID, data = pet.save, FUN = function(x) max(x, na.rm = TRUE))
petSD <- aggregate(. ~ ID, data = pet.save, FUN = function(x) sd(x, na.rm = TRUE))

PET <- cbind(PETavg = PETavg[,2], PETmin = PETmin[,2], PETmax = PETmax[,2], petSD = petSD[,2])
rownames(PET) <- sp.names


##########
# elev
#########
elev.vals <- extract(elev, distEA, na.rm=T)
str(elev.vals)

elev.save <- elev.vals

#mean by groups - ID 
ELEVavg <- aggregate(. ~ ID, data = elev.save, FUN = function(x) mean(x, na.rm = TRUE))
ELEVmin <- aggregate(. ~ ID, data = elev.save, FUN = function(x) min(x, na.rm = TRUE))
ELEVmax <- aggregate(. ~ ID, data = elev.save, FUN = function(x) max(x, na.rm = TRUE))
elevSD <- aggregate(. ~ ID, data = elev.save, FUN = function(x) sd(x, na.rm = TRUE))

ELEV <- cbind(ELEVavg = ELEVavg[,2], ELEVmin = ELEVmin[,2], ELEVmax = ELEVmax[,2], elevSD = elevSD[,2])
rownames(ELEV) <- sp.names


#########
# tmax
#########
tmax.vals <- extract(tmax, distEA, na.rm=T)
str(tmax.vals)

tmax.save <- tmax.vals

#mean by groups - ID 
result <- aggregate(. ~ ID, data = tmax.save, FUN = function(x) mean(x, na.rm = TRUE))
tmax.mean <- result

unique(result$ID)

sp.names <- values(distEA)
sp.names <- sp.names$sci_name

rownames(tmax.mean) <- sp.names
tmax.mean <- tmax.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
TMAXavg <- apply(tmax.mean, 1, mean, na.rm=T)
tmaxSD <- apply(tmax.mean, 1, sd, na.rm=T)

tmax_total <- data.frame(TMAXavg, tmaxSD)


#########
# tmin
#########
tmin.vals <- extract(tmin, distEA, na.rm=T)
str(tmin.vals)

tmin.save <- tmin.vals

#mean by groups - ID 
result <- aggregate(. ~ ID, data = tmin.save, FUN = function(x) mean(x, na.rm = TRUE))
tmin.mean <- result

rownames(tmin.mean) <- sp.names
tmin.mean <- tmin.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
TMINavg <- apply(tmin.mean, 1, mean, na.rm=T)
tminSD <- apply(tmin.mean, 1, sd, na.rm=T)

tmin_total <- data.frame(TMINavg, tminSD)


########
# TAVG
########
tavg.vals <- extract(tavg, distEA, na.rm=T)
str(tavg.vals)

tavg.save <- tavg.vals

#mean by groups - ID 
result <- aggregate(. ~ ID, data = tavg.save, FUN = function(x) mean(x, na.rm = TRUE))
tavg.mean <- result

rownames(tavg.mean) <- sp.names
tavg.mean <- tavg.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
TAVGavg <- apply(tavg.mean, 1, mean, na.rm=T)
tavgSD <- apply(tavg.mean, 1, sd, na.rm=T)

tavg_total <- data.frame(TAVGavg, tavgSD)


########
# VAPR
########
vapr.vals <- extract(vapr, distEA, na.rm=T)
str(vapr.vals)

vapr.save <- vapr.vals

result <- aggregate(. ~ ID, data = vapr.save, FUN = function(x) mean(x, na.rm = TRUE))
vapr.mean <- result

rownames(vapr.mean) <- sp.names
vapr.mean <- vapr.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
VAPRavg <- apply(vapr.mean, 1, mean, na.rm=T)
VAPRmax <- apply(vapr.mean, 1, max, na.rm=T)
VAPRmin <- apply(vapr.mean, 1, min, na.rm=T)
vaprSD <- apply(vapr.mean, 1, sd, na.rm=T)

vapr_total <- data.frame(VAPRavg, vaprSD)


###########
# WIND
##########
wind.vals <- extract(wind, distEA, na.rm=T)
str(wind.vals)

wind.save <- wind.vals

result <- aggregate(. ~ ID, data = wind.save, FUN = function(x) mean(x, na.rm = TRUE))
wind.mean <- result

rownames(wind.mean) <- sp.names
wind.mean <- wind.mean[-1]

# With that we can calculate maximum and minimum "mean" NDVI
WINDavg <- apply(wind.mean, 1, mean, na.rm=T)
WINDmax <- apply(wind.mean, 1, max, na.rm=T)
WINDmin <- apply(wind.mean, 1, min, na.rm=T)
windSD <- apply(wind.mean, 1, sd, na.rm=T)

wind_total <- data.frame(WINDavg, windSD)


# IMPORTANT
# NOTE: for precipitation and solar radiation we will take a different approach:
# first we will add all the layers to estimate the anual precipitation and s. radiation

###########
# PREC
##########
# first sum all the layers
sum.prec <- sum(prec)

prec.vals <- extract(sum.prec, distEA, na.rm=T)
str(prec.vals)

prec.save <- prec.vals

PRECavg <- aggregate(. ~ ID, data = prec.save, FUN = function(x) mean(x, na.rm = TRUE))
PRECmin <- aggregate(. ~ ID, data = prec.save, FUN = function(x) min(x, na.rm = TRUE))
PRECmax <- aggregate(. ~ ID, data = prec.save, FUN = function(x) max(x, na.rm = TRUE))
precSD <- aggregate(. ~ ID, data = prec.save, FUN = function(x) sd(x, na.rm = TRUE))

PREC <- cbind(PRECavg = PRECavg[,2], PRECmin = PRECmin[,2], PRECmax = PRECmax[,2], precSD = precSD[,2])
rownames(PREC) <- sp.names



#####################
# srad
#####################
sum.srad <- sum(srad)

srad.vals <- extract(sum.srad, distEA, na.rm=T)
str(srad.vals)

srad.save <- srad.vals

SRADavg <- aggregate(. ~ ID, data = srad.save, FUN = function(x) mean(x, na.rm = TRUE))
SRADmin <- aggregate(. ~ ID, data = srad.save, FUN = function(x) min(x, na.rm = TRUE))
SRADmax <- aggregate(. ~ ID, data = srad.save, FUN = function(x) max(x, na.rm = TRUE))
sradSD <- aggregate(. ~ ID, data = srad.save, FUN = function(x) sd(x, na.rm = TRUE))


##########################################
# Organize and save the data
##########################################

env95 <- as.data.frame(cbind(ARIDavg[,2], ARIDmax[,2], ARIDmin[,2], aridSD[,2],
                             NDVIavg, NDVImax, NDVImin, ndviSD,
                             PETavg[,2], PETmax[,2], PETmin[,2], petSD[,2],
                             PRECavg[,2], PRECmax[,2], PRECmin[,2], precSD[,2], 
                             SRADavg[,2], SRADmax[,2], SRADmin[,2], sradSD[,2],
                             TAVGavg, tavgSD, 
                             TMAXavg, tmaxSD, 
                             TMINavg, tminSD, 
                             VAPRavg, VAPRmax, VAPRmin, vaprSD,
                             ELEVavg[,2], ELEVmax[,2], ELEVmin[,2], elevSD[,2],
                             WINDavg, WINDmax, WINDmin, windSD))

colnames(env95) <- c("ARIDavg", "ARIDmax", "ARIDmin", "aridSD", "NDVIavg", "NDVImax", "NDVImin", "ndviSD",
                     "PETavg",  "PETmax", "PETmin", "petSD", "PRECavg", "PRECmax", "PRECmin", "precSD",
                     "SRADavg", "SRADmax", "SRADmin", "sradSD", "TAVGavg", "tavgSD", "TMAXavg", "tmaxSD",
                     "TMINavg", "tminSD", "VAPRavg", "VAPRmax", "VAPRmin", "vaprSD", "ELEVavg", "ELEVmax",
                     "ELEVmin", "elevSD", "WINDavg", "WINDmax", "WINDmin", "windSD")

rownames(env95) <- distEA$sci_name
rownames(env95)

# write.csv(env95, "env_vars_EA_sp.csv")


