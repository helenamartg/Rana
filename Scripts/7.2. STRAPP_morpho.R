#################################################
# STRAPP on the rest of morhpological variables
#################################################
# For this analysis we need the output of BAMM

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

library(phytools)
library(ape)
library(BAMMtools)


# BAMM output diversification
tr59 <- read.tree("Data/phyloRana_59sp.tre")
ediv <- getEventData(tr59, eventdata = "BAMM/event_data_div59.txt", 
                     burnin=0.1, type = "diversification")


#1. HL 
#=====
tr.hl<-read.tree("BAMM/tree_hl.tre")
edata_hl<- getEventData(tr.hl, eventdata = "BAMM/event_data_trait_hl.txt", 
                        burnin=0.5, type = "trait")


# rates of phenotypic evolution and diversification 
hl <- getTipRates(edata_hl) 
div <- getTipRates(ediv, returnNetDiv = T)

hl.rates <- hl$beta.avg  # Average tip phenotypic rates
div.rates <- div$netdiv.avg      # Average tip diversification rates

# remove sp
names(hl.rates) == names(div.rates)
sp <- as.factor(intersect(names(hl.rates), names(div.rates)))
div.rates <- as.data.frame(div.rates)
div.rates$sp <- rownames(div.rates)

div.rates.hl <- div.rates[names(hl.rates),]
div.rates.hl <- setNames(div.rates.hl$div.rates, div.rates.hl$sp)
names(hl.rates) == names(div.rates.hl)

# Extract from diversification analyses the sp for which we have data of HL
subphy <- subtreeBAMM(ediv,tips=names(hl.rates))
subphy$tip.label == names(hl.rates)
hl.rates <- as.data.frame(hl.rates)
hl.rates$sp <- rownames(hl.rates)
hl.rates <- hl.rates[subphy$tip.label,]
hl.rates$sp==subphy$tip.label
hl.rates <- setNames(hl.rates$hl.rates, hl.rates$sp)
subphy$tip.label == names(hl.rates)

hl.permu.sp<-traitDependentBAMM(subphy, hl.rates, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

hl.permu.sp$estimate # Very low correlation between speciation rates and svl evolution
hl.permu.sp$p.value  # NOT significant


# Separate radiations
par(mar=c(1,1,1,1))
plot(tr.hl)
nodelabels(cex=0.5)

# AMERICAN
#=========
am <- extract.clade(tr.hl, node=56)
am$tip.label

hl.rates <-as.data.frame(hl.rates)
hl.rates$species <- rownames(hl.rates)
hl.rates.am <- hl.rates[am$tip.label,]

subphy.am <- subtreeBAMM(ediv, tips=hl.rates.am$species)
subphy.am$tip.label == hl.rates.am$species
hl.rates.am <- hl.rates.am[subphy.am$tip.label,]
hl.rates.am$sp==subphy.am$tip.label

hl.rates.am <- setNames(hl.rates.am$hl.rates, hl.rates.am$species)
subphy.am$tip.label == names(hl.rates.am)


hl.permu.am<-traitDependentBAMM(subphy.am, hl.rates.am, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

hl.permu.am$estimate # Very low correlation between speciation rates and svl evolution
hl.permu.am$p.value  # NOT significant


# EURASIAN
#==========
plot(tr.hl)
nodelabels(cex=0.5)
ea <- drop.clade(tr.hl, am$tip.label)
ea <- drop.tip(ea, "NA")
plot(ea)

ea$tip.label
hl.rates.ea <- hl.rates[ea$tip.label,]

subphy.ea <- subtreeBAMM(ediv, tips=hl.rates.ea$species)
subphy.ea$tip.label == hl.rates.ea$species
hl.rates.ea <- hl.rates.ea[subphy.ea$tip.label,]
hl.rates.ea$sp==subphy.ea$tip.label

hl.rates.ea <- setNames(hl.rates.ea$hl.rates, hl.rates.ea$species)
subphy.ea$tip.label == names(hl.rates.ea)


hl.permu.ea<-traitDependentBAMM(subphy.ea, hl.rates.ea, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

hl.permu.ea$estimate # Very low correlation between speciation rates and svl evolution
hl.permu.ea$p.value  # NOT significant

cor(hl.rates.ea, getTipRates(subphy.ea, returnNetDiv = T)$netdiv.avg)
anova(lm(getTipRates(subphy.ea, returnNetDiv = T)$netdiv.avg ~ hl.rates.ea))

#PLOT 3 HL (GGPLOT)
div.ea <- getTipRates(subphy.ea, returnNetDiv = T)
div.ea.rates <- div.ea$netdiv.avg

div.am <- getTipRates(subphy.am, returnNetDiv = T)
div.am.rates <- div.am$netdiv.avg

library(ggplot2)
mat.all.hl <- as.data.frame(cbind(hl.rates, div.rates.hl))
mat.ea.hl <- as.data.frame(cbind(hl.rates.ea, div.ea.rates))
mat.am.hl <- as.data.frame(cbind(hl.rates.am, div.am.rates))

p1 <- ggplot(mat.all.hl, aes(x=hl.rates, y=div.rates.hl)) + 
  geom_point() +
  theme_classic() +
  ylab("Diversification Rate") +
  xlab("Rate of HL evolution") +
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Global")

p2 <- ggplot(mat.ea.hl, aes(x=hl.rates.ea, y=div.ea.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of HL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Eurasia")

p3 <- ggplot(mat.am.hl, aes(x=hl.rates.am, y=div.am.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of HL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") +
  ggtitle("America")

library(patchwork)
p1 + p2 + p3



#########
# 2. HW 
#########
tr.hw<-read.tree("BAMM/tree_hw.tre")
edata_hw<- getEventData(tr.hw, eventdata = "BAMM/event_data_trait_hw.txt", 
                        burnin=0.5, type = "trait")


# rates of phenotypic evolution and diversification 
hw <- getTipRates(edata_hw) 
hw.rates <- hw$beta.avg  # Average tip phenotypic rates


# remove sp
names(hw.rates) == names(div.rates)
sp <- as.factor(intersect(names(hw.rates), names(div.rates)))
div.rates <- as.data.frame(div.rates)
div.rates$sp <- rownames(div.rates)

div.rates.hw <- div.rates[names(hw.rates),]
div.rates.hw <- setNames(div.rates.hw$div.rates, div.rates.hw$sp)
names(hw.rates) == names(div.rates.hw)



# Extract from diversification analyses the sp for which we have data of HL
subphy <- subtreeBAMM(ediv,tips=names(hw.rates))
subphy$tip.label == names(hw.rates)
hw.rates <- as.data.frame(hw.rates)
hw.rates$sp <- rownames(hw.rates)
hw.rates <- hw.rates[subphy$tip.label,]
hw.rates$sp==subphy$tip.label
hw.rates <- setNames(hw.rates$hw.rates, hw.rates$sp)
subphy$tip.label == names(hw.rates)

hw.permu.sp<-traitDependentBAMM(subphy, hw.rates, reps=10000, return.full = T,  
                                two.tailed = T, method='s')


hw.permu.sp$estimate # Very low correlation between speciation rates and svl evolution
hw.permu.sp$p.value  # NOT significant


# Separate radiations
par(mar=c(1,1,1,1))
plot(tr.hw)
nodelabels(cex=0.5)

# AMERICAN
#=========
am <- extract.clade(tr.hw, node=44)
am$tip.label

hw.rates <-as.data.frame(hw.rates)
hw.rates$species <- rownames(hw.rates)
hw.rates.am <- hw.rates[am$tip.label,]

subphy.am <- subtreeBAMM(ediv, tips=hw.rates.am$species)
subphy.am$tip.label == hw.rates.am$species

hw.rates.am <- setNames(hw.rates.am$hw.rates, hw.rates.am$species)
subphy.am$tip.label == names(hw.rates.am)

hw.permu.am<-traitDependentBAMM(subphy.am, hw.rates.am, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

hw.permu.am$estimate # Very low correlation between speciation rates and svl evolution
hw.permu.am$p.value  # NOT significant


# EURASIAN
#==========
plot(tr.hw)
nodelabels(cex=0.5)
ea <- drop.clade(tr.hw, am$tip.label)
ea <- drop.tip(ea, "NA")
plot(ea)

ea$tip.label
hw.rates.ea <- hw.rates[ea$tip.label,]

subphy.ea <- subtreeBAMM(ediv, tips=hw.rates.ea$species)
subphy.ea$tip.label == hw.rates.ea$species
# hw.rates.ea <- hw.rates.ea[subphy.ea$tip.label,]
# hw.rates.ea$sp==subphy.ea$tip.label

hw.rates.ea <- setNames(hw.rates.ea$hw.rates, hw.rates.ea$species)
subphy.ea$tip.label == names(hw.rates.ea)

hw.permu.ea<-traitDependentBAMM(subphy.ea, hw.rates.ea, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

hw.permu.ea$estimate # Very low correlation between speciation rates and svl evolution
hw.permu.ea$p.value  # NOT significant


#PLOT 3 HL (GGPLOT)
div.ea <- getTipRates(subphy.ea, returnNetDiv = T)
div.ea.rates <- div.ea$netdiv.avg

div.am <- getTipRates(subphy.am, returnNetDiv = T)
div.am.rates <- div.am$netdiv.avg

library(ggplot2)
mat.all.hw <- as.data.frame(cbind(hw.rates, div.rates.hw))
mat.ea.hw <- as.data.frame(cbind(hw.rates.ea, div.ea.rates))
mat.am.hw <- as.data.frame(cbind(hw.rates.am, div.am.rates))

p4 <- ggplot(mat.all.hw, aes(x=hw.rates, y=div.rates.hw)) + 
  geom_point() +
  theme_classic() +
  ylab("Diversification Rate") +
  xlab("Rate of HW evolution") +
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Global")

p5 <- ggplot(mat.ea.hw, aes(x=hw.rates.ea, y=div.ea.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of HW evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Eurasia")

p6 <- ggplot(mat.am.hw, aes(x=hw.rates.am, y=div.am.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of HW evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") +
  ggtitle("America")

library(patchwork)
p4 + p5 + p6




#########
# 3. FL 
#########
tr.fl<-read.tree("BAMM/tree_fl.tre")
edata_fl<- getEventData(tr.fl, eventdata = "BAMM/event_data_trait_fl.txt", 
                        burnin=0.5, type = "trait")


# rates of phenotypic evolution and diversification 
fl <- getTipRates(edata_fl) 
fl.rates <- fl$beta.avg  # Average tip phenotypic rates

# remove sp
names(fl.rates) == names(div.rates)
sp <- as.factor(intersect(names(fl.rates), names(div.rates)))
# div.rates <- as.data.frame(div.rates)
# div.rates$sp <- rownames(div.rates)

div.rates.fl <- div.rates[names(fl.rates),]
div.rates.fl <- setNames(div.rates.fl$div.rates, div.rates.fl$sp)
names(fl.rates) == names(div.rates.fl)

# Extract from diversification analyses the sp for which we have data of HL
subphy <- subtreeBAMM(ediv,tips=names(fl.rates))
subphy$tip.label == names(fl.rates)
fl.rates <- as.data.frame(fl.rates)
fl.rates$sp <- rownames(fl.rates)
fl.rates <- fl.rates[subphy$tip.label,]
fl.rates$sp==subphy$tip.label
fl.rates <- setNames(fl.rates$fl.rates, fl.rates$sp)
subphy$tip.label == names(fl.rates)

fl.permu.sp<-traitDependentBAMM(subphy, fl.rates, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

fl.permu.sp$estimate # Very low correlation between speciation rates and svl evolution
fl.permu.sp$p.value  # NOT significant


# Separate radiations
plot(tr.fl)
nodelabels(cex=0.5)

# AMERICAN
#=========
am <- extract.clade(tr.fl, node=53)
am$tip.label

fl.rates <-as.data.frame(fl.rates)
fl.rates$species <- rownames(fl.rates)
fl.rates.am <- fl.rates[am$tip.label,]

subphy.am <- subtreeBAMM(ediv, tips=fl.rates.am$species)
subphy.am$tip.label == fl.rates.am$species
fl.rates.am <- fl.rates.am[subphy.am$tip.label,]
fl.rates.am$sp==subphy.am$tip.label

fl.rates.am <- setNames(fl.rates.am$fl.rates, fl.rates.am$species)
subphy.am$tip.label == names(fl.rates.am)

fl.permu.am<-traitDependentBAMM(subphy.am, fl.rates.am, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

fl.permu.am$estimate # Very low correlation between speciation rates and svl evolution
fl.permu.am$p.value  # NOT significant


# EURASIAN
#==========
ea <- drop.clade(tr.fl, am$tip.label)
ea <- drop.tip(ea, "NA")
plot(ea)

ea$tip.label
fl.rates.ea <- fl.rates[ea$tip.label,]

subphy.ea <- subtreeBAMM(ediv, tips=fl.rates.ea$species)
subphy.ea$tip.label == fl.rates.ea$species
fl.rates.ea <- fl.rates.ea[subphy.ea$tip.label,]
fl.rates.ea$sp==subphy.ea$tip.label

fl.rates.ea <- setNames(fl.rates.ea$fl.rates, fl.rates.ea$species)
subphy.ea$tip.label == names(fl.rates.ea)

fl.permu.ea<-traitDependentBAMM(subphy.ea, fl.rates.ea, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

fl.permu.ea$estimate # Very low correlation between speciation rates and svl evolution
fl.permu.ea$p.value  # NOT significant


#PLOT 3 FL (GGPLOT)
div.ea <- getTipRates(subphy.ea, returnNetDiv = T)
div.ea.rates <- div.ea$netdiv.avg

div.am <- getTipRates(subphy.am, returnNetDiv = T)
div.am.rates <- div.am$netdiv.avg

library(ggplot2)
mat.all.fl <- as.data.frame(cbind(fl.rates, div.rates.fl))
mat.ea.fl <- as.data.frame(cbind(fl.rates.ea, div.ea.rates))
mat.am.fl <- as.data.frame(cbind(fl.rates.am, div.am.rates))

p7 <- ggplot(mat.all.fl, aes(x=fl.rates, y=div.rates.fl)) + 
  geom_point() +
  theme_classic() +
  ylab("Diversification Rate") +
  xlab("Rate of fL evolution") +
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Global")

p8 <- ggplot(mat.ea.fl, aes(x=fl.rates.ea, y=div.ea.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of FL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Eurasia")

p9 <- ggplot(mat.am.fl, aes(x=fl.rates.am, y=div.am.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of FL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") +
  ggtitle("America")

library(patchwork)
p7 + p8 + p9


#########
# 4. TL 
#########
tr.tl<-read.tree("BAMM/tree_tl.tre")
edata_tl<- getEventData(tr.tl, eventdata = "BAMM/event_data_trait_tl.txt", 
                        burnin=0.5, type = "trait")

tl <- getTipRates(edata_tl) 
tl.rates <- tl$beta.avg  # Average tip phenotypic rates

# remove sp
names(tl.rates) == names(div.rates)
sp <- as.factor(intersect(names(tl.rates), names(div.rates)))

div.rates.tl <- div.rates[names(tl.rates),]
div.rates.tl <- setNames(div.rates.tl$div.rates, div.rates.tl$sp)
names(tl.rates) == names(div.rates.tl)

# Extract from diversification analyses the sp for which we have data of HL
subphy <- subtreeBAMM(ediv,tips=names(tl.rates))
subphy$tip.label == names(tl.rates)
tl.rates <- as.data.frame(tl.rates)
tl.rates$sp <- rownames(tl.rates)
tl.rates <- tl.rates[subphy$tip.label,]
tl.rates$sp==subphy$tip.label
tl.rates <- setNames(tl.rates$tl.rates, tl.rates$sp)
subphy$tip.label == names(tl.rates)

tl.permu.sp<-traitDependentBAMM(subphy, tl.rates, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

tl.permu.sp$estimate # Very low correlation between speciation rates and svl evolution
tl.permu.sp$p.value  # NOT significant


# Separate radiations
plot(tr.tl)
nodelabels(cex=0.5)

# AMERICAN
#=========
am <- extract.clade(tr.tl, node=57)
am$tip.label

tl.rates <-as.data.frame(tl.rates)
tl.rates$species <- rownames(tl.rates)
tl.rates.am <- tl.rates[am$tip.label,]

subphy.am <- subtreeBAMM(ediv, tips=tl.rates.am$species)
subphy.am$tip.label == tl.rates.am$species
tl.rates.am <- tl.rates.am[subphy.am$tip.label,]
tl.rates.am$sp==subphy.am$tip.label

tl.rates.am <- setNames(tl.rates.am$tl.rates, tl.rates.am$species)
subphy.am$tip.label == names(tl.rates.am)

tl.permu.am<-traitDependentBAMM(subphy.am, tl.rates.am, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

tl.permu.am$estimate # Very low correlation between speciation rates and svl evolution
tl.permu.am$p.value  # NOT significant


# EURASIAN
#==========
ea <- drop.clade(tr.tl, am$tip.label)
ea <- drop.tip(ea, "NA")
plot(ea)

ea$tip.label
tl.rates.ea <- tl.rates[ea$tip.label,]

subphy.ea <- subtreeBAMM(ediv, tips=tl.rates.ea$species)
subphy.ea$tip.label == tl.rates.ea$species
tl.rates.ea <- tl.rates.ea[subphy.ea$tip.label,]
tl.rates.ea$sp==subphy.ea$tip.label

tl.rates.ea <- setNames(tl.rates.ea$tl.rates, tl.rates.ea$species)
subphy.ea$tip.label == names(tl.rates.ea)

tl.permu.ea<-traitDependentBAMM(subphy.ea, tl.rates.ea, reps=10000, return.full = T,  
                                two.tailed = T, method='s')

tl.permu.ea$estimate # Very low correlation between speciation rates and svl evolution
tl.permu.ea$p.value  # NOT significant


#PLOT 3 TL (GGPLOT)
div.ea <- getTipRates(subphy.ea, returnNetDiv = T)
div.ea.rates <- div.ea$netdiv.avg

div.am <- getTipRates(subphy.am, returnNetDiv = T)
div.am.rates <- div.am$netdiv.avg

library(ggplot2)
mat.all.tl <- as.data.frame(cbind(tl.rates, div.rates.tl))
mat.ea.tl <- as.data.frame(cbind(tl.rates.ea, div.ea.rates))
mat.am.tl <- as.data.frame(cbind(tl.rates.am, div.am.rates))

p10 <- ggplot(mat.all.tl, aes(x=tl.rates, y=div.rates.tl)) + 
  geom_point() +
  theme_classic() +
  ylab("Diversification Rate") +
  xlab("Rate of TL evolution") +
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Global")

p11 <- ggplot(mat.ea.tl, aes(x=tl.rates.ea, y=div.ea.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of TL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") + 
  ggtitle("Eurasia")

p12 <- ggplot(mat.am.tl, aes(x=tl.rates.am, y=div.am.rates)) + 
  geom_point() +
  theme_classic() +
  ylab("")+
  xlab("Rate of TL evolution")+
  geom_smooth(method=lm, se=FALSE, color="red") +
  ggtitle("America")

library(patchwork)
p10 + p11 + p12



# MEGAPLOTTTT
library(ggpubr)
ggarrange(pa, pb, pc,
          p1, p2, p3, p4, p5, p6,
          p7, p8, p9, p10, p11, p12,
          ncol = 3, nrow = 5)

