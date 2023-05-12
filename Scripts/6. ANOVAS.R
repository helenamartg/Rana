###################
# ANOVAS of SVL
#################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

library(RRPP)
library(ape)
library(phytools)

# 1. Import data
#================
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
dt74 <- dt74[,-1]

dt68 <-read.csv("DATA/dt_68sp.csv")
dt68 <- dt68[,-1]

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
dtEA <- dtEA[match(dtEA$sp, phyloEA$tip.label),]
dtAM <- dtAM[match(dtAM$sp, phyloAM$tip.label),]
dt74 <- dt74[match(dt74$sp, phyloRana$tip.label),]


# 5. Group by clutch size (new categorical variable)
#====================================================
range(na.omit(dt74$clutch_size))

# Small: 100 - 999
# Medium: 1000 - 4999
# Large: 5000 - 17000

dt74$clutch <- NA
for (i in 1:nrow(dt74)){
  if (is.na(dt74[i,7])){
    dt74$clutch[i] <- NA
  } else if (dt74[i,7] < 999){
    dt74$clutch[i] <- "Small"
    } else if (dt74[i,7] > 999 & dt74[i,7] <= 4999){
      dt74$clutch[i] <- "Medium"
     } else if (dt74[i,7] > 4999){
        dt74$clutch[i] <- "Large"
     }
}

table(dt74$clutch)


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


rrpp.ea <- rrpp.data.frame(svl = log(dtEA$SVL), 
                           sp = as.factor(dtEA$sp),
                           gen = as.factor(dtEA$genus),
                           ecotype = dtEA$ecotype,
                           state = dtEA$state,
                           radiation = dtEA$phylo_reg,
                           biome = dtEA$biome,
                           clutch = log(dtEA$clutch_size),
                           phy.cvc = cov.ea)

rrpp.am <- rrpp.data.frame(svl = log(dtAM$SVL), 
                           sp = as.factor(dtAM$sp),
                           gen = as.factor(dtAM$genus),
                           ecotype = dtAM$ecotype,
                           state = dtAM$state,
                           radiation = dtAM$phylo_reg,
                           biome = dtAM$biome,
                           clutch = log(dtAM$clutch_size),
                           phy.cvc = cov.am)

rrpp.global <- rrpp.data.frame(svl = log(dt74$SVL), 
                               sp = as.factor(dt74$sp),
                               gen = as.factor(dt74$genus),
                               ecotype = dt74$ecotype,
                               state = dt74$state,
                               radiation = dt74$phylo_reg,
                               biome = dt74$biome,
                               clutch = log(dt74$clutch_size),
                               phy.cvc = cov.rana)


# Are body sizes different across radiations?
lm.SVL <- lm.rrpp(svl ~ radiation, data = rrpp.global, SS.type = "III")  # Without phylogeny (GLS)
anova(lm.SVL)   # Significant

lm.SVLphylo <- lm.rrpp(svl ~ radiation, data = rrpp.global, SS.type = "III", Cov = cov.rana)  # phylogeny (PGLS)
anova(lm.SVLphylo)   # Not significant


par(mar=c(4,4,4,4))
boxplot(rrpp.global$svl ~ rrpp.global$radiation, col=c("gold", "darkorchid"), ylab = "log(SVL)", xlab="Radiation", names=c("American", "Eurasian"))
stripchart(rrpp.global$svl ~ rrpp.global$radiation, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")


# Are body sizes different across simple states?
lm.state <- lm.rrpp(svl ~ state, data = rrpp.global, SS.type = "III") # Without phylogeny (GLS)
anova(lm.state) # Significant

lm.statephylo <- lm.rrpp(svl ~ state, data = rrpp.global, SS.type = "III", Cov = cov.rana) #  phylogeny (PGLS)
anova(lm.statephylo) # Significant

pal <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
boxplot(rrpp.global$svl ~ rrpp.global$state, col=pal, xlab="", ylab="log(SVL)")
stripchart(rrpp.global$svl ~ rrpp.global$state, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")


# Are body sizes different across ecotypes?
lm.eco <- lm.rrpp(svl ~ ecotype, data = rrpp.global, SS.type = "III") # Non Phylogenetic (GLS)
anova(lm.eco)

lm.ecophylo <- lm.rrpp(svl ~ ecotype, data = rrpp.global, SS.type = "III", Cov = cov.rana) # Phylogenetic (PGLS)
anova(lm.ecophylo)

boxplot(rrpp.global$svl ~ rrpp.global$ecotype, col=c("dodgerblue4", "cornflowerblue", "chocolate4"))
stripchart(rrpp.global$svl ~ rrpp.global$ecotype, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")


# Are body sizes different across biomes?
lm.biome <- lm.rrpp(svl ~ biome, data=rrpp.global, SS.type="III") # Non Phylogenetic (GLS)
anova(lm.biome)

lm.biomephylo <- lm.rrpp(svl ~ biome, data=rrpp.global, SS.type="III", Cov = cov.rana) # Phylogenetic (PGLS)
anova(lm.biomephylo) # Significant

par(mar=c(4,4,4,4))
pal2 <- c("cadetblue3", "antiquewhite3", "chocolate", "darkolivegreen4", "chartreuse3")
boxplot(rrpp.global$svl ~ rrpp.global$biome, xlab="", ylab="log(SVL)", col=pal2)
stripchart(rrpp.global$svl ~ rrpp.global$biome, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")

table(rrpp.global$biome)

# Is body size related with clutch size
lm.clutch <- lm.rrpp(svl ~ clutch, data=rrpp.global, Cov = cov.rana, SS.type = "III")
anova(lm.clutch) # Significant
plot(rrpp.global$clutch, rrpp.global$svl, pch=21, bg="black", 
     ylab="svl", xlab="clutch", main="All species")
abline(lm.clutch, lwd=2, col="red")


# INTERACTIONS
# Effect of radiation in the relationship between svl and region
lm.radstate <- lm.rrpp(svl ~ radiation*state, data = rrpp.global, SS.type = "III", Cov = cov.rana)
anova(lm.radstate)

# Effect of the ecotype in svl ~ biome
lm.ecob <- lm.rrpp(svl ~ ecotype*biome, data = rrpp.global, SS.type = "III", Cov = cov.rana)
anova(lm.ecob)

# Effect of radiation in svl ~ ecotype
lm.ecorad <- lm.rrpp(svl ~ radiation*ecotype, data = rrpp.global, SS.type = "III", Cov = cov.rana)
anova(lm.ecorad)

# Effect of radiation in svl ~ biome
lm.radb <- lm.rrpp(svl ~ radiation*biome, data = rrpp.global, SS.type = "III", Cov = cov.rana)
anova(lm.radb)

table(interaction(rrpp.global$ecotype, rrpp.global$biome))
table(rrpp.global$ecotype)
table(interaction(rrpp.global$biome, rrpp.global$radiation))



# 7. SEPARATED RADIATIONS
#==========================
# Here we take into account phylogeny (PGLS)
##############
# 7.1 Eurasia
##############

# Is SVL different across biomes?
lm.EA1 <- lm.rrpp(svl ~ biome, data = rrpp.ea, SS.type = "III", Cov = cov.ea)
anova(lm.EA1)  #No significant differences

boxplot(rrpp.ea$svl ~ rrpp.ea$biome, xlab="", ylab="log(SVL)", 
        col=c("cadetblue3", "chocolate", "darkolivegreen4", "chartreuse3"), main="Eurasia")
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
anova(lm.EA3) #Significant differences

boxplot(rrpp.ea$svl ~ rrpp.ea$state, xlab="", ylab="log(SVL)", 
        col=c("darkorchid1", "darkorchid4"), main="Eurasia")
stripchart(rrpp.ea$svl ~ rrpp.ea$state, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")

# Is body size related with clutch size?
lm.EA4 <- lm.rrpp(svl ~ clutch, data=rrpp.ea, Cov = cov.ea, SS.type = "III")
anova(lm.EA4) #significant

plot(rrpp.ea$clutch, rrpp.ea$svl, pch=21, bg="black", 
     main="European radiation", ylab="log(SVL)", xlab="log(Clutch size)")
abline(lm.EA4, col="darkorchid", lwd=2)


##############
# 7.2 America
##############
lm.AM1 <- lm.rrpp(svl ~ biome, data = rrpp.am, SS.type = "III", Cov = cov.am)
anova(lm.AM1)  #Significant differences

boxplot(rrpp.am$svl ~ rrpp.am$biome, xlab="", ylab="log(SVL)", col=pal2, main="America")
stripchart(rrpp.am$svl ~ rrpp.am$biome, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")

lm.AM2 <- lm.rrpp(svl ~ ecotype, data = rrpp.am, SS.type = "III", Cov=cov.am)
anova(lm.AM2)  # No significant differences

boxplot(rrpp.am$svl ~ rrpp.am$ecotype, xlab="", ylab="log(SVL)", 
        col=c("dodgerblue4", "cornflowerblue", "chocolate4"), main="America")
stripchart(rrpp.am$svl ~ rrpp.am$ecotype, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")

lm.AM3 <- lm.rrpp(svl ~ state, data=rrpp.am, SS.type = "III", Cov = cov.am)
anova(lm.AM3) #Not sign differences

boxplot(rrpp.am$svl ~ rrpp.am$state, xlab="", ylab="log(SVL)", 
        col=c("gold4", "gold3", "gold"), main="America")
stripchart(rrpp.am$svl ~ rrpp.am$state, vertical = T, method = "jitter", pch = 19, 
           add = TRUE, col = "black")

lm.AM4 <- lm.rrpp(svl ~clutch, data=rrpp.am, Cov = cov.am, SS.type = "III")
anova(lm.AM4)

plot(rrpp.am$clutch, rrpp.am$svl, pch=21, bg="black", 
     main="American radiation", ylab="log(SVL)", xlab="log(Clutch size)")
abline(lm.AM4, col="gold", lwd=2)
