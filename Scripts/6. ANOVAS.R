###################
# ANOVAS of SVL
#################
# 74 species for body size 
# 68 for morphological variables
rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(RRPP)
library(ape)
library(phytools)
library(castor)

# 1. Import data
#================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
dt74 <- dt74[,-1]

dt68 <-read.csv("DATA/dt_68sp.csv")
dt68 <- dt68[,-1]
rownames(dt68) <- dt68$sp


# 2. Import phylogeny
#=====================
phyloRana <- read.tree("DATA/phyloRana_74sp.tre")
phyloRana$tip.label
phyloRana_68 <- read.tree("DATA/phyloRana_68sp.tre")
phyloEA <- read.tree("DATA/phyloEA_31sp.tre")
phyloAM <- read.tree("DATA/phyloAM_37sp.tre")


# 3. Split dataset
#==================
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]

# 4. Match with phylogeny
#=========================
dtEA <- dtEA[phyloEA$tip.label,]
dtAM <- dtAM[phyloAM$tip.label,]
dt74 <- dt74[phyloRana$tip.label,]
dt68 <-  dt68[phyloRana_68$tip.label,]

dt74$sp == phyloRana$tip.label
dt68$sp == phyloRana_68$tip.label

dtEA$sp == phyloEA$tip.label
dtAM$sp == phyloAM$tip.label



# 6. GLOBAL ANOVAS with SVL
#==========================
# We need to create three datasets: 
# 1) one with 31 sp from eurasian radiation --> rrpp.ea
# 2) one with 37sp from american radiation --> rrpp.am
# 3) another one with 74 sp for global comparisions --> rrpp.global

# Build the Phylogenetic var-covariance matrix (correlation)
cov.ea <- vcv(phyloEA, model="Brownian") 
cov.am <- vcv(phyloAM, model="Brownian") 
cov.rana <- vcv(phyloRana, model="Brownian")
cov.rana68 <- vcv(phyloRana_68, model="Brownian")


rrpp.ea <- rrpp.data.frame(svl = dtEA$log_svl, 
                           sp = as.factor(dtEA$sp),
                           ecotype = as.factor(dtEA$ecotype),
                           state = as.factor(dtEA$state),
                           biome = as.factor(dtEA$biome),
                           phy.cvc = cov.ea)

rrpp.am <- rrpp.data.frame(svl = dtAM$log_svl, 
                           sp = as.factor(dtAM$sp),
                           ecotype = as.factor(dtAM$ecotype),
                           state = as.factor(dtAM$state),
                           biome = as.factor(dtAM$biome),
                           phy.cvc = cov.am)

rrpp.global <- rrpp.data.frame(svl = dt74$log_svl, 
                               sp = as.factor(dt74$sp),
                               ecotype = as.factor(dt74$ecotype),
                               state = as.factor(dt74$state),
                               radiation = as.factor(dt74$phylo_reg),
                               biome = as.factor(dt74$biome),
                               phy.cvc = cov.rana)

rrpp.global68 <- rrpp.data.frame(svl = dt68$log_svl, 
                               sp = as.factor(dt68$sp),
                               ecotype = as.factor(dt68$ecotype),
                               state = as.factor(dt68$state),
                               radiation = as.factor(dt68$phylo_reg),
                               biome = as.factor(dt68$biome),
                               phy.cvc = cov.rana68)


# Are body sizes different across radiations?
lm.SVL <- lm.rrpp(svl ~ radiation, data = rrpp.global, SS.type = "III")  # Without phylogeny (GLS)
anova(lm.SVL)   # Significant

lm.SVLphylo <- lm.rrpp(svl ~ radiation, data = rrpp.global, SS.type = "III", Cov = cov.rana)  # phylogeny (PGLS)
anova(lm.SVLphylo)   # Not significant

lm.SVLphylo68 <- lm.rrpp(svl ~ radiation, data=rrpp.global68, SS.type = "III", Cov = cov.rana68)
anova(lm.SVLphylo68)  # Not significant

par(mar=c(4,4,4,4), cex.axis=1.3)
boxplot(rrpp.global68$svl ~ rrpp.global68$radiation, col=c("gold", "darkorchid"), ylab = "log(SVL)", xlab="Radiation", names=c("American", "Eurasian"), cex.lab=1.5)
stripchart(rrpp.global$svl ~ rrpp.global$radiation, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black", cex=0.5)



# Are body sizes different across simple states?
lm.state <- lm.rrpp(svl ~ state, data = rrpp.global, SS.type = "III") # Without phylogeny (GLS)
anova(lm.state) # Significant

lm.statephylo <- lm.rrpp(svl ~ state, data = rrpp.global, SS.type = "III", Cov = cov.rana) #PGLS
anova(lm.statephylo) # Not Significant

lm.statephylo68 <- lm.rrpp(svl ~ state, data = rrpp.global68, SS.type = "III", Cov = cov.rana68)
anova(lm.statephylo68) # Not signif

pal <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
boxplot(rrpp.global$svl ~ rrpp.global$state, col=pal, xlab="", ylab="log(SVL)")
stripchart(rrpp.global$svl ~ rrpp.global$state, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black", cex=0.5)



# Are body sizes different across ecotypes?
lm.eco <- lm.rrpp(svl ~ ecotype, data = rrpp.global, SS.type = "III") # Non Phylogenetic (GLS)
anova(lm.eco)

lm.ecophylo <- lm.rrpp(svl ~ ecotype, data = rrpp.global, SS.type = "III", Cov = cov.rana) # PGLS
anova(lm.ecophylo) # significant

lm.ecophylo68 <- lm.rrpp(svl ~ ecotype, data = rrpp.global68, SS.type = "III", Cov = cov.rana68) # PGLS
anova(lm.ecophylo68) # Not significant

boxplot(rrpp.global$svl ~ rrpp.global$ecotype, col=c("dodgerblue4", "cornflowerblue", "chocolate4"))
stripchart(rrpp.global$svl ~ rrpp.global$ecotype, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black", cex=0.5)




# Are body sizes different across biomes?
lm.biome <- lm.rrpp(svl ~ biome, data=rrpp.global, SS.type="III") # Non Phylogenetic (GLS)
anova(lm.biome)

lm.biomephylo <- lm.rrpp(svl ~ biome, data=rrpp.global, SS.type="III", Cov = cov.rana) # PGLS
anova(lm.biomephylo) # Significant

lm.biomephylo68 <- lm.rrpp(svl ~ biome, data=rrpp.global68, SS.type="III", Cov = cov.rana68) # PGLS
anova(lm.biomephylo68) # Significant


par(mar=c(4,4,4,4), cex.axis=0.8)
pal2 <- c("cadetblue3", "antiquewhite3", "chocolate", "darkolivegreen4", "chartreuse3")
boxplot(rrpp.global$svl ~ rrpp.global$biome, xlab="", ylab="log(SVL)", col=pal2,
        names=c("Boreal", "Desert", "Grassland", "Temperate", "Tropical"))
stripchart(rrpp.global$svl ~ rrpp.global$biome, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black", cex=0.5)


# Is body size related with clutch size?
# lm.clutch <- lm.rrpp(svl ~ clutch, data=rrpp.global, Cov = cov.rana, SS.type = "III")
# anova(lm.clutch) # Significant
# plot(rrpp.global$clutch, rrpp.global$svl, pch=21, bg="black", 
#      ylab="svl", xlab="clutch", main="All species")
# abline(lm.clutch, lwd=2, col="red")

# Create a table for SAVE
gls.svl <- rbind(anova(lm.SVL)$table, anova(lm.state)$table, anova(lm.biome)$table, anova(lm.eco)$table)
pgls.svl <- rbind(anova(lm.SVLphylo)$table, anova(lm.statephylo)$table, anova(lm.biomephylo)$table,
                  anova(lm.ecophylo)$table)
pgls.svl.68 <- rbind(anova(lm.SVLphylo68)$table, anova(lm.statephylo68)$table, anova(lm.biomephylo68)$table,
                     anova(lm.ecophylo68)$table)
anovas.svl <- cbind(gls.svl, "", pgls.svl, "", pgls.svl.68)
# write.table(anovas.svl, "Results/Anovas_svl.csv")



# 7. SEPARATED RADIATIONS
#==========================
# Here we take into account phylogeny (PGLS)

##################
# 7.1 SVL Eurasia
##################
# Is SVL different across biomes?
lm.EA1 <- lm.rrpp(svl ~ biome, data = rrpp.ea, SS.type = "III", Cov = cov.ea)
anova(lm.EA1)  #No significant differences

par(cex.axis=0.8)
boxplot(rrpp.ea$svl ~ rrpp.ea$biome, xlab="", ylab="log(SVL)", 
        col=c("cadetblue3", "chocolate", "darkolivegreen4", "chartreuse3"), main="Eurasia",
        names=c("Boreal", "Grassland", "Temperate", "Tropical"))
stripchart(rrpp.ea$svl ~ rrpp.ea$biome, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")


# Is SVL different across ecotypes?
lm.EA2 <- lm.rrpp(svl ~ ecotype, data = rrpp.ea, SS.type = "III",  Cov = cov.ea)
anova(lm.EA2) # No significant differences

boxplot(rrpp.ea$svl ~ rrpp.ea$ecotype, xlab="", ylab="log(SVL)", 
        col=c("dodgerblue4", "cornflowerblue", "chocolate4"), main="Eurasia")
stripchart(rrpp.ea$svl ~ rrpp.ea$ecotype, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")


# Is SVL different between states?
lm.EA3 <- lm.rrpp(svl ~ state, data = rrpp.ea, SS.type = "III",  Cov = cov.ea)
anova(lm.EA3) #No significant differences

boxplot(rrpp.ea$svl ~ rrpp.ea$state, xlab="", ylab="log(SVL)", 
        col=c("darkorchid1", "darkorchid4"), main="Eurasia")
points(tapply(rrpp.ea$svl, rrpp.ea$state, mean), pch=4, lwd=2)

svl.ea <- rbind(anova(lm.EA3)$table, anova(lm.EA2)$table, anova(lm.EA1)$table)


##################
# 7.2 SVL America
##################
# Biome
# Remove R_sylvatica because is the only one in Boreal forest
dtAM <- dtAM[dtAM$sp!="Rana_sylvatica",]
phyloAM2 <- get_subtree_with_tips(phyloRana, dtAM$sp)$subtree
cov.am <- vcv(phyloAM2, model="Brownian")
rrpp.am <- rrpp.data.frame(svl = dtAM$log_svl, 
                           sp = as.factor(dtAM$sp),
                           ecotype = as.factor(dtAM$ecotype),
                           state = as.factor(dtAM$state),
                           biome = as.factor(dtAM$biome),
                           phy.cvc = cov.am)

lm.AM1 <- lm.rrpp(svl ~ biome, data = rrpp.am, SS.type = "III", Cov = cov.am)
anova(lm.AM1)  #Significant differences

par(cex.axis=0.8)
boxplot(rrpp.am$svl ~ rrpp.am$biome, xlab="", ylab="log(SVL)", col=pal2[-1], main="America",
        names=c("Desert","Grassland", "Temperate", "Tropical"))
points(tapply(rrpp.am$svl, rrpp.am$biome, mean), pch=4, lwd=2)

# Ecotype
lm.AM2 <- lm.rrpp(svl ~ ecotype, data = rrpp.am, SS.type = "III", Cov=cov.am)
anova(lm.AM2)  # No significant differences

boxplot(rrpp.am$svl ~ rrpp.am$ecotype, xlab="", ylab="log(SVL)", 
        col=c("dodgerblue4", "cornflowerblue", "chocolate4"), main="America")


# Simple state
lm.AM3 <- lm.rrpp(svl ~ state, data=rrpp.am, SS.type = "III", Cov = cov.am)
anova(lm.AM3) #Not sign differences

boxplot(rrpp.am$svl ~ rrpp.am$state, xlab="", ylab="log(SVL)", 
        col=c("gold4", "gold3", "gold"), main="America")

svl.am <- rbind(anova(lm.AM3)$table, anova(lm.AM2)$table, anova(lm.AM1)$table)

# SAVE
# write.table(svl.am, "Results/Anovas_svl_AM.csv")
# write.table(svl.ea, "Results/Anovas_svl_EA.csv")



# 8. MORPHOLOGY
#===============
# Remaining morphological variables
morpho <- read.csv("DATA/dt_morpho_55sp.csv")
rownames(morpho) <- morpho$sp

#remove 6 species from Eurasia
morpho <- morpho[morpho$sp!="Rana_boylii",]
morpho <- morpho[morpho$sp!="Rana_sierrae",]
morpho <- morpho[morpho$sp!="Rana_muscosa",]


# Add column with biomes
biomes <- read.csv("DATA/raw_data/biome_all.csv")
rownames(biomes) <- biomes$sp
biomes <- as.data.frame(biomes)
biomes$sp <- rownames(biomes)
biomes <- biomes[,5:6]

biom <- NULL
for (i in 1:nrow(biomes)){
  for (j in 1:nrow(morpho)){
    if (biomes[i,2]== morpho[j,2]){
      x <- biomes[i, c(1:2)]
      biom <- rbind(biom, x)
    }
  }
}
dim(biom)
dim(morpho)
morpho$sp==biom$sp
biom <- biom[morpho$sp,]

morpho <- cbind(morpho, biom)
morpho <- morpho[morpho$sp!= "Rana_tlaloci",]

morpho.res <- morpho[, c(2,4:6, 13:17,19)]



##############################
#      REDUCED DATASET 
# (remove all rows with NAs)
##############################
morpho.res
red.morpho <- na.omit(morpho.res) #34sp
phylo34 <- get_subtree_with_tips(phyloRana, red.morpho$sp)$subtree

# Check order
red.morpho <- red.morpho[phylo34$tip.label,]
red.morpho$sp == phylo34$tip.label

# Cov matrix
cov.rana34 <- vcv(phylo34, model = "Brownian")

# Loop PGLS
mat <- NULL
for (i in 6:9){
  temp <- rrpp.data.frame(morpho = red.morpho[,i],
                          state = red.morpho$state,
                          radiation = red.morpho$phylo_reg,
                          biome = red.morpho$new_biome,
                          ecotype = red.morpho$ecotype,
                          Cov = cov.rana34)
  lm.rad <- lm.rrpp(morpho ~ radiation, data = temp, SS.type = "III", Cov = cov.rana34)  
  anova(lm.rad)

  lm.state <- lm.rrpp(morpho ~ state, data = temp, SS.type = "III", Cov = cov.rana34)
  anova(lm.state)

  lm.ecotype <- lm.rrpp(morpho ~ ecotype, data = temp, SS.type = "III", Cov = cov.rana34)
    anova(lm.ecotype)

  lm.biome <- lm.rrpp(morpho~ biome, data = temp, SS.type = "III", Cov = cov.rana34)
  anova(lm.biome)
  
  lms <- rbind(anova(lm.rad)$table, anova(lm.state)$table,
                 anova(lm.ecotype)$table, anova(lm.biome)$table)

  lms <- cbind(lms, names(red.morpho[i]))

  mat <- rbind(mat, lms)
}

anovas_morpho_all <- mat


# EURASIA
#==========
morpho.ea <- red.morpho[red.morpho$phylo_reg=="EA",] # 21 sp
phylo.ea <- get_subtree_with_tips(phyloEA, morpho.ea$sp)$subtree

morpho.ea$sp == phylo.ea$tip.label

cov.ea <- vcv(phylo.ea, model = "Brownian")

mat <- NULL
for (i in 6:9){
  temp <- rrpp.data.frame(morpho = morpho.ea[,i],
                          state = morpho.ea$state,
                          biome = morpho.ea$new_biome,
                          ecotype = morpho.ea$ecotype,
                          Cov = cov.ea)
  
  lm.state <- lm.rrpp(morpho ~ state, data = temp, SS.type = "III", Cov = cov.ea)
  anova(lm.state)
  
  lm.ecotype <- lm.rrpp(morpho ~ ecotype, data = temp, SS.type = "III", Cov = cov.ea)
  anova(lm.ecotype)
  
  lm.biome <- lm.rrpp(morpho~ biome, data = temp, SS.type = "III", Cov = cov.ea)
  anova(lm.biome)
  
  lms <- rbind(anova(lm.state)$table,
               anova(lm.ecotype)$table, anova(lm.biome)$table)
  
  lms <- cbind(lms, names(morpho.ea[i]))
  
  mat <- rbind(mat, lms)
}

anovas_morpho_EA <- mat



# AMERICA
#==========
morpho.am <- red.morpho[red.morpho$phylo_reg=="AM",] # 13 sp
phylo.am <- get_subtree_with_tips(phyloAM, morpho.am$sp)$subtree

morpho.am$sp == phylo.am$tip.label

cov.am <- vcv(phylo.am, model = "Brownian")

mat <- NULL
for (i in 6:9){
  temp <- rrpp.data.frame(morpho = morpho.am[,i],
                          state = morpho.am$state,
                          biome = morpho.am$new_biome,
                          ecotype = morpho.am$ecotype,
                          Cov = cov.am)
  
  lm.state <- lm.rrpp(morpho ~ state, data = temp, SS.type = "III", Cov = cov.am)
  anova(lm.state)
  
  lm.ecotype <- lm.rrpp(morpho ~ ecotype, data = temp, SS.type = "III", Cov = cov.am)
  anova(lm.ecotype)
  
  lm.biome <- lm.rrpp(morpho~ biome, data = temp, SS.type = "III", Cov = cov.am)
  anova(lm.biome)
  
  lms <- rbind(anova(lm.state)$table,
               anova(lm.ecotype)$table, anova(lm.biome)$table)
  
  lms <- cbind(lms, names(morpho.am[i]))
  
  mat <- rbind(mat, lms)
}

anovas_morpho_AM <- mat


# SAVE
# write.table(anovas_morpho_all, "Results/ANOVAS/Anovas_morpho_all.csv")
# write.table(anovas_morpho_EA, "Results/ANOVAS/Anovas_morpho_EA.csv")
# write.table(anovas_morpho_AM, "Results/ANOVAS/Anovas_morpho_AM.csv")







# A PARTIR DE AQUI FATUL################################


# MEGALOOP to not miss species and clip dataset and phylogeny by each trait
mat <- NULL
for (i in 6:9){
  morpho.x <- morpho.res[complete.cases(morpho.res[,i]),]      # Remove sp from morpho dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.x$sp))$subtree  # Clip phylogeny

  morpho.x <- morpho.x[phylo.x$tip.label,]

  cov.temp <- vcv(phylo.x, model="Brownian")

  rrpp.morpho <- rrpp.data.frame(morpho = morpho.x[,i],
                  sp = as.factor(morpho.x$sp),
                  ecotype = as.factor(morpho.x$ecotype),
                  state = as.factor(morpho.x$state),
                  radiation = as.factor(morpho.x$phylo_reg),
                  biome = as.factor(morpho.x$new_biome),
                  phy.cvc = cov.temp)

  lm.rad <- lm.rrpp(morpho ~ radiation, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.rad)

  lm.state <- lm.rrpp(morpho ~ state, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.state)

  lm.ecotype <- lm.rrpp(morpho ~ ecotype, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.ecotype)

  lm.biome <- lm.rrpp(morpho~ biome, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.biome)

 #MERGE
  names(morpho.x)
  lms <- rbind(anova(lm.rad)$table, anova(lm.state)$table, anova(lm.biome)$table)

  lms <- cbind(lms, names(morpho.x[i]))

  mat <- rbind(mat, lms)

}

# SAVE
# write.table(mat, "Results/Anovas_morpho_clippingfilo.csv")


# 8.1 Morpho EURASIA
######################
morpho.res
# Remove Amerana subgenera
morpho.res <- morpho.res[morpho.res$sp!="Rana_boylii",]
morpho.res <- morpho.res[morpho.res$sp!="Rana_sierrae",]
morpho.res <- morpho.res[morpho.res$sp!="Rana_muscosa",]

morpho.ea <- morpho.res[morpho.res$phylo_reg=="EA",]

mat.ea <- NULL
for (i in 6:9){
  morpho.x <- morpho.ea[complete.cases(morpho.ea[,i]),]      # Remove sp from morpho dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.x$sp))$subtree  # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  
  cov.temp <- vcv(phylo.x, model="Brownian") 
  
  rrpp.morpho <- rrpp.data.frame(morpho = morpho.x[,i], 
                                 sp = as.factor(morpho.x$sp),
                                 ecotype = as.factor(morpho.x$ecotype),
                                 state = as.factor(morpho.x$state),
                                 radiation = as.factor(morpho.x$phylo_reg),
                                 biome = as.factor(morpho.x$new_biome),
                                 phy.cvc = cov.temp)
  
  lm.state <- lm.rrpp(morpho ~ state, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.state) 
  
  lm.ecotype <- lm.rrpp(morpho ~ ecotype, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.ecotype)
  
  lm.biome <- lm.rrpp(morpho~ biome, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.biome) 
  
  #MERGE
  names(morpho.x)
  lms <- rbind(anova(lm.state)$table, 
               anova(lm.ecotype)$table,
               anova(lm.biome)$table)
  
  lms <- cbind(lms, names(morpho.x[i]))
  
  mat.ea <- rbind(mat.ea, lms)
  
}

# SAVE
# write.table(mat.ea, "Results/ANOVAS/Anovas_morpho_EA_clippingfilo.csv")


# 8.2 Morpho AMERICA
######################
morpho.am <- morpho.res[morpho.res$phylo_reg=="AM",]

mat.am <- NULL
for (i in 6:9){
  morpho.x <- morpho.am[complete.cases(morpho.am[,i]),]      # Remove sp from morpho dataset
  phylo.x <- (get_subtree_with_tips(phyloRana, morpho.x$sp))$subtree  # Clip phylogeny
  
  morpho.x <- morpho.x[phylo.x$tip.label,]
  
  cov.temp <- vcv(phylo.x, model="Brownian") 
  
  rrpp.morpho <- rrpp.data.frame(morpho = morpho.x[,i], 
                                 sp = as.factor(morpho.x$sp),
                                 ecotype = as.factor(morpho.x$ecotype),
                                 state = as.factor(morpho.x$state),
                                 radiation = as.factor(morpho.x$phylo_reg),
                                 biome = as.factor(morpho.x$new_biome),
                                 phy.cvc = cov.temp)
  
  lm.state <- lm.rrpp(morpho ~ state, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.state) 
  
  lm.ecotype <- lm.rrpp(morpho ~ ecotype, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.ecotype)
  
  lm.biome <- lm.rrpp(morpho~ biome, data = rrpp.morpho, SS.type = "III", Cov = cov.temp)
  anova(lm.biome) 
  
  #MERGE
  names(morpho.x)
  lms <- rbind(anova(lm.state)$table, 
               anova(lm.ecotype)$table, anova(lm.biome)$table)
  
  lms <- cbind(lms, names(morpho.x[i]))
  
  mat.am <- rbind(mat.am, lms)
  
}

# SAVE
# write.table(mat.am, "Results/ANOVAS/Anovas_morpho_AM_clippingfilo.csv")



# 9. MANOVAS   (PABLO)
#=======================
morpho2 <- na.omit(morpho.res) 
tree2 <- get_subtree_with_tips(phyloRana, rownames(morpho2))$subtree
state <- setNames(as.factor(morpho2$state), morpho2$sp)
biomes <- setNames(as.factor(morpho2$new_biome), morpho2$sp)
ecotypes <- setNames(as.factor(morpho2$ecotype), morpho2$sp)
morpho <- as.data.frame(morpho2[, 6:9])
radiation <- setNames(as.factor(morpho2$phylo_reg), morpho2$sp)

class(state)

# PCA
pca.glob <- prcomp(morpho, scale. = T, center = T)
summary(pca.glob)
cor(morpho, pca.glob$x)  

pc <- as.data.frame(pca.glob$x[,1:2])


#Is general morphology different between states?
phyMANOVA.1 <- geiger::aov.phylo(morpho ~ state, tree2, nsim=50, test="Wilks")   # significant ¿?
summary.manova(phyMANOVA.1, test="Wilks")  # lambda wilks = 0.0726 P-value significant 

phyMANOVA.11 <- geiger::aov.phylo(pc ~ state, tree2, nsim=50, test="Wilks") # signf considering phylogeny
summary(phyMANOVA.11, test="Wilks") # Significant


#And between biomes?
phyMANOVA.2 <- geiger::aov.phylo(morpho ~ biomes, tree2, nsim=50, test="Wilks")  # Significant
summary(phyMANOVA.2)

phyMANOVA.22 <- geiger::aov.phylo(pc ~ biomes, tree2, nsim=50, test="Wilks") #Not signf considering phylogeny
summary(phyMANOVA.22, test="Wilks") # Significant 

# And between radiaions?
phyMANOVA.3 <- geiger::aov.phylo(morpho ~ radiation, tree2, nsim=50, test="Wilks")  #  NOT Significant
summary.manova(phyMANOVA.3)

phyMANOVA.33 <- geiger::aov.phylo(pc ~ radiation, tree2, nsim=50, test="Wilks") # signf considering phylogeny
summary(phyMANOVA.33, test="Wilks") # Significant

# And ecotypes?
phyMANOVA.4 <- geiger::aov.phylo(morpho ~ ecotypes, tree2, nsim=50, test="Wilks")  #  Non Significant
summary.manova(phyMANOVA.4)

phyMANOVA.44 <- geiger::aov.phylo(pc ~ ecotypes, tree2, nsim=50, test="Wilks") # signf considering phylogeny
summary(phyMANOVA.44, test="Wilks") # Significant



# EUROPE
morpho_EA <- morpho2[morpho2$phylo_reg=="EA",]
tr_EA <- get_subtree_with_tips(phyloRana, rownames(morpho_EA))$subtree
state_EA <- setNames(as.factor(morpho_EA$state), morpho_EA$sp)
biomes_EA <- setNames(as.factor(morpho_EA$new_biome), morpho_EA$sp)
ecotypes_EA <- setNames(as.factor(morpho_EA$ecotype), morpho_EA$sp)
morpho_EA <- as.data.frame(morpho_EA[, 6:9])


phyMANOVA.EA <- geiger::aov.phylo(morpho_EA ~ state_EA, tr_EA, nsim=50, test="Wilks")   # significant ¿?
summary.manova(phyMANOVA.EA, test="Wilks")  # lambda wilks = 0.0726 P-value significant 

phyMANOVA.EA <- geiger::aov.phylo(morpho_EA ~ biomes_EA, tr_EA, nsim=50, test="Wilks")   # significant ¿?
summary.manova(phyMANOVA.EA, test="Wilks") 


#AMERICAN
morpho_AM <- morpho2[morpho2$phylo_reg=="AM",]
tr_AM <- get_subtree_with_tips(phyloRana, rownames(morpho_AM))$subtree
state_AM <- setNames(as.factor(morpho_AM$state), morpho_AM$sp)
biomes_AM <- setNames(as.factor(morpho_AM$new_biome), morpho_AM$sp)
ecotypes_AM <- setNames(as.factor(morpho_AM$ecotype), morpho_AM$sp)
morpho_AM <- as.data.frame(morpho_AM[, 6:9])

phyMANOVA.AM <- geiger::aov.phylo(morpho_AM ~ state_AM, tr_AM, nsim=50, test="Wilks")   # significant ¿?
summary.manova(phyMANOVA.AM, test="Wilks")

phyMANOVA.AM <- geiger::aov.phylo(morpho_AM ~ biomes_AM, tr_AM, nsim=50, test="Wilks")   # significant ¿?
summary.manova(phyMANOVA.AM, test="Wilks")
