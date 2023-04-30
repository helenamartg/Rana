####################################
#         IUCN ANURAN 
# Script to obtain distributions
####################################


rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/My analyses/1. Distributions")

library(maptools)
library(gdata)
library(rworldmap)
library(raster)
library(rgdal)

amphibians <- readShapeSpatial('ANURA/ANURA.shp', delete_null_obj=T)

# get all the Rana spp?
Ranidae.spp <- amphibians[amphibians@data$family=='RANIDAE',]
Ranidae.spp@data <- drop.levels(Ranidae.spp@data)
genus <- unique(substr(Ranidae.spp@data$sci_name, 1, regexpr(' ', Ranidae.spp@data$sci_name)-1))


# Extract genera
Rana.spp <- Ranidae.spp[Ranidae.spp@data$genus=='Rana',]
Rana.spp@data <- drop.levels(Rana.spp@data)
Rana.spp@data$sci_name

Lithobates.spp <- amphibians[amphibians@data$genus=='Lithobates',]
Lithobates.spp@data <- drop.levels(Lithobates.spp@data)
Lithobates.spp@data$sci_name

Pseudo.spp <- amphibians[amphibians@data$genus=='Pseudorana',]
Pseudo.spp@data <- drop.levels(Pseudo.spp@data)
Pseudo.spp@data$sci_name

# plot all species' distributions by genus
# newmap <- getMap(resolution = "low")


# Rana in phylo (37 spp)
sort(levels(Rana.spp@data$sci_name))
# 1,2,3,4,5,6,7,8,10,11,13,14,15,17,18,20,21,22,23,24,25,26,
# 27,29,30,31,32,35,36,37,38,39,40,42,43,44,45
Rsp <- sort(levels(Rana.spp@data$sci_name))
subRana <- c(1,2,3,4,5,6,7,8,10,13,15,16,17,19,20,22,24,26,
  27,28,30,31,32,34,35,37,40,42,43,44,45,46,47,49,50,52)
Rspp <- Rsp[subRana]

# Lithobates in phylo (39 spp)
sort(levels(Lithobates.spp@data$sci_name))
# 1,2,3,4,5,6,8,9,10,12,13,14,16,18,19,20,23,24,25,26,27,28,
# 29,30,31,33,34,35,36,37,38,40,41,43,44,45,47,48,49
Lsp <- sort(levels(Lithobates.spp@data$sci_name))
subLitho <- c(1,2,3,5,6,7,9,10,11,13,14,15,17,21,22,23,26,27,28,
  29,30,31,32,33,34,36,37,38,39,40,41,42,43,45,46,47,49,50,51)
Lspp <- Lsp[subLitho]

# Pseudorana in phylo (1 sp)
Psp <- sort(levels(Pseudo.spp@data$sci_name))
# 2 (i.e. Pseudorana weiningensis)


Rana.spp.phylo <- Rana.spp[Rana.spp@data$sci_name==Rspp[1] |
    Rana.spp@data$sci_name==Rspp[2] | Rana.spp@data$sci_name==Rspp[3] |
    Rana.spp@data$sci_name==Rspp[4] | Rana.spp@data$sci_name==Rspp[5] |
    Rana.spp@data$sci_name==Rspp[6] | Rana.spp@data$sci_name==Rspp[7] |
    Rana.spp@data$sci_name==Rspp[8] | Rana.spp@data$sci_name==Rspp[9] |
    Rana.spp@data$sci_name==Rspp[10] | Rana.spp@data$sci_name==Rspp[11] |
    Rana.spp@data$sci_name==Rspp[12] | Rana.spp@data$sci_name==Rspp[13] |
    Rana.spp@data$sci_name==Rspp[14] | Rana.spp@data$sci_name==Rspp[15] |
    Rana.spp@data$sci_name==Rspp[16] | Rana.spp@data$sci_name==Rspp[17] |
    Rana.spp@data$sci_name==Rspp[18] | Rana.spp@data$sci_name==Rspp[19] |
    Rana.spp@data$sci_name==Rspp[20] | Rana.spp@data$sci_name==Rspp[21] |
    Rana.spp@data$sci_name==Rspp[22] | Rana.spp@data$sci_name==Rspp[23] |
    Rana.spp@data$sci_name==Rspp[24] | Rana.spp@data$sci_name==Rspp[25] |
    Rana.spp@data$sci_name==Rspp[26] | Rana.spp@data$sci_name==Rspp[27] |
    Rana.spp@data$sci_name==Rspp[28] | Rana.spp@data$sci_name==Rspp[29] |
    Rana.spp@data$sci_name==Rspp[30] | Rana.spp@data$sci_name==Rspp[31] |
    Rana.spp@data$sci_name==Rspp[32] | Rana.spp@data$sci_name==Rspp[33] |
    Rana.spp@data$sci_name==Rspp[34] | Rana.spp@data$sci_name==Rspp[35] |
    Rana.spp@data$sci_name==Rspp[36],]

Rana.spp.phylo@data <- drop.levels(Rana.spp.phylo@data)
unique(Rana.spp.phylo@data$sci_name)



Litho.spp.phylo <- Lithobates.spp[Lithobates.spp@data$sci_name==Lspp[1] |
    Lithobates.spp@data$sci_name==Lspp[2] | Lithobates.spp@data$sci_name==Lspp[3] |
    Lithobates.spp@data$sci_name==Lspp[4] | Lithobates.spp@data$sci_name==Lspp[5] |
    Lithobates.spp@data$sci_name==Lspp[6] | Lithobates.spp@data$sci_name==Lspp[7] |
    Lithobates.spp@data$sci_name==Lspp[8] | Lithobates.spp@data$sci_name==Lspp[9] |
    Lithobates.spp@data$sci_name==Lspp[10] | Lithobates.spp@data$sci_name==Lspp[11] |
    Lithobates.spp@data$sci_name==Lspp[12] | Lithobates.spp@data$sci_name==Lspp[13] |
    Lithobates.spp@data$sci_name==Lspp[14] | Lithobates.spp@data$sci_name==Lspp[15] |
    Lithobates.spp@data$sci_name==Lspp[16] | Lithobates.spp@data$sci_name==Lspp[17] |
    Lithobates.spp@data$sci_name==Lspp[18] | Lithobates.spp@data$sci_name==Lspp[19] |
    Lithobates.spp@data$sci_name==Lspp[20] | Lithobates.spp@data$sci_name==Lspp[21] |
    Lithobates.spp@data$sci_name==Lspp[22] | Lithobates.spp@data$sci_name==Lspp[23] |
    Lithobates.spp@data$sci_name==Lspp[24] | Lithobates.spp@data$sci_name==Lspp[25] |
    Lithobates.spp@data$sci_name==Lspp[26] | Lithobates.spp@data$sci_name==Lspp[27] |
    Lithobates.spp@data$sci_name==Lspp[28] | Lithobates.spp@data$sci_name==Lspp[29] |
    Lithobates.spp@data$sci_name==Lspp[30] | Lithobates.spp@data$sci_name==Lspp[31] |
    Lithobates.spp@data$sci_name==Lspp[32] | Lithobates.spp@data$sci_name==Lspp[33] |
    Lithobates.spp@data$sci_name==Lspp[34] | Lithobates.spp@data$sci_name==Lspp[35] |
    Lithobates.spp@data$sci_name==Lspp[36] | Lithobates.spp@data$sci_name==Lspp[37] |
    Lithobates.spp@data$sci_name==Lspp[38] | Lithobates.spp@data$sci_name==Lspp[39],]

Litho.spp.phylo@data <- drop.levels(Litho.spp.phylo@data)
unique(Litho.spp.phylo@data$sci_name)


Pseudo.spp.phylo <- Pseudo.spp[Pseudo.spp@data$sci_name=='Pseudorana weiningensis',]
Pseudo.spp.phylo@data <- drop.levels(Pseudo.spp.phylo@data)
unique(Pseudo.spp.phylo@data$sci_name)


distRana <- rbind(Rana.spp.phylo, Litho.spp.phylo, Pseudo.spp.phylo)
sort(unique(distRana@data$sci_name))

# Pseudorana 
distEA <- rbind(Rana.spp.phylo, Pseudo.spp.phylo)
distAM <- Litho.spp.phylo



#FIRST REMOVE INTRODUCED POPS
levels(distRana@data$legend)
distRana <- distRana[distRana@data$origin==1,] # remove introduced populations
distRana <- distRana[distRana@data$presence == 1,] # remove extinct, presence uncertain, and everything that is not Exant (resident)

distRana@data$presence

# change names to match phylogeny
unique(distRana@data$sci_name)
distRana@data$sci_name <- as.character(distRana@data$sci_name)

for(i in 1:length(distRana@data$sci_name)){
  
  if(substr(distRana@data$sci_name[i], 1, regexpr(' ',distRana@data$sci_name[i])-1)=='Rana'){
    distRana@data$sci_name[i] <- gsub(' ','_',distRana@data$sci_name[i])
  }
  
  if(substr(distRana@data$sci_name[i], 1, regexpr(' ',distRana@data$sci_name[i])-1)=='Lithobates'){
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='areolatus'){
      distRana@data$sci_name[i] <- 'Rana_areolata'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='neovolcanicus'){
      distRana@data$sci_name[i] <- 'Rana_neovolcanica'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='omiltemanus'){
      distRana@data$sci_name[i] <- 'Rana_omiltemana'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='maculatus'){
      distRana@data$sci_name[i] <- 'Rana_maculata'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='pustulosus'){
      distRana@data$sci_name[i] <- 'Rana_pustulosa'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='sevosus'){
      distRana@data$sci_name[i] <- 'Rana_sevosa'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='catesbeianus'){
      distRana@data$sci_name[i] <- 'Rana_catesbeiana'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='sphenocephalus'){
      distRana@data$sci_name[i] <- 'Rana_sphenocephala'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='sylvaticus'){
      distRana@data$sci_name[i] <- 'Rana_sylvatica'
    }
    
    if(substr(distRana@data$sci_name[i], regexpr(' ',distRana@data$sci_name[i])+1, nchar(distRana@data$sci_name[i][1]))=='vibicarius'){
      distRana@data$sci_name[i] <- 'Rana_vibicaria'
    }
    
    else{
      distRana@data$sci_name[i] <- gsub('Lithobates ','Rana_',distRana@data$sci_name[i])
    }
    
  }
  
  if(substr(distRana@data$sci_name[i], 1, regexpr(' ',distRana@data$sci_name[i])-1)=='Pseudorana'){
    distRana@data$sci_name[i] <- gsub('Pseudorana ','Rana_',distRana@data$sci_name[i])
  }
}
  
distRana@data$sci_name <- as.factor(distRana@data$sci_name)
sort(unique(distRana@data$sci_name))


# same thing for the European radiation
distEA <- distEA[distEA@data$origin==1,] # remove introduced populations
distEA <- distEA[distEA@data$presence == 1,]

unique(distEA@data$sci_name)
distEA@data$sci_name <- as.character(distEA@data$sci_name)

for(i in 1:length(distEA@data$sci_name)){
  
  distEA@data$sci_name[i] <- gsub(' ','_',distEA@data$sci_name[i])
}


distEA@data$sci_name <- as.factor(distEA@data$sci_name)
sort(unique(distEA@data$sci_name))



# same thing for the American radiation
distAM <- distAM[distAM@data$origin==1,] # remove introduced populations
distAM <- distAM[distAM@data$presence == 1,]

unique(distAM@data$sci_name)
distAM@data$sci_name <- as.character(distAM@data$sci_name)

for(i in 1:length(distAM@data$sci_name)){
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='areolatus'){
    distAM@data$sci_name[i] <- 'Rana_areolata'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='neovolcanicus'){
    distAM@data$sci_name[i] <- 'Rana_neovolcanica'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='omiltemanus'){
    distAM@data$sci_name[i] <- 'Rana_omiltemana'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='maculatus'){
    distAM@data$sci_name[i] <- 'Rana_maculata'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='pustulosus'){
    distAM@data$sci_name[i] <- 'Rana_pustulosa'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='sevosus'){
    distAM@data$sci_name[i] <- 'Rana_sevosa'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='catesbeianus'){
    distAM@data$sci_name[i] <- 'Rana_catesbeiana'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='sphenocephalus'){
    distAM@data$sci_name[i] <- 'Rana_sphenocephala'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='sylvaticus'){
    distAM@data$sci_name[i] <- 'Rana_sylvatica'
  }
  
  if(substr(distAM@data$sci_name[i], regexpr(' ',distAM@data$sci_name[i])+1, nchar(distAM@data$sci_name[i][1]))=='vibicarius'){
    distAM@data$sci_name[i] <- 'Rana_vibicaria'
  }
  
  else{
    distAM@data$sci_name[i] <- gsub('Lithobates ','Rana_',distAM@data$sci_name[i])
  }
}


distAM@data$sci_name <- as.factor(distAM@data$sci_name)
sort(unique(distAM@data$sci_name))


# Clean the environment
# rm(list=keep(distRana, distEA, distAM))

# CHECK
plot(distEA[distEA@data$sci_name=='Rana_iberica',]) # still good!

###### Now we want to have one single polygon per species
# EURASIA
library(rgeos)
distEA2 <- gUnaryUnion(distEA, id = distEA@data$sci_name, checkValidity = T)
distEA2 <- as(distEA2, "SpatialPolygonsDataFrame")
distEA2@data$sci_name <- row.names(distEA2)

# more checking
plot(distEA2[distEA2@data$sci_name=='Rana_sierrae',]) # still good!
plot(distEA2[distEA2@data$sci_name=='Rana_iberica',])
plot(distEA2[distEA2@data$sci_name=='Rana_temporaria',])
plot(distEA2[distEA2@data$sci_name=='Rana_amurensis',])

# SAVE
# fold <- paste(getwd(), '/dist_Eurasia', sep='')
# writeOGR(distEA2, fold, "dist_Eurasia", driver="ESRI Shapefile", overwrite=T)


# AMERICA
distAM2 <- gUnaryUnion(distAM, id = distAM@data$id_no)
distAM2 <- as(distAM2, "SpatialPolygonsDataFrame")
distAM2@data$sci_name <- row.names(distAM2)

plot(distAM2[distAM2@data$sci_name=='Rana_catesbeiana',]) # still good!
plot(distAM2[distAM2@data$sci_name=='Rana_areolata',]) # still good!

# SAVE
# fold <- paste(getwd(), '/dist_America', sep='')
# writeOGR(distAM2, fold, "dist_America", driver="ESRI Shapefile", overwrite=T)


