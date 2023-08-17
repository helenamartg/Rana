##########################
# GGPLOTS Anovas
##########################

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

# 2. Split dataset
#==================
dtEA <- dt68[dt68$phylo_reg=="EA",]
dtAM <- dt68[dt68$phylo_reg=="AM",]


# 3. GGPLOT: SVL
#================
library(ggplot2)
library(ggtree)

#3.1 Global
############
pal <- c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")
pal2 <- c("cadetblue3", "antiquewhite3", "chocolate", "darkolivegreen4", "chartreuse3")

svl.rad <- ggplot(dt74, aes(y=log_svl, x=reorder(real_reg, log_svl))) + 
  stat_boxplot(geom = "errorbar",width = 0.25) +
  geom_boxplot(color="black", fill=c("darkorchid", "gold3"), alpha=1, varwidth = T)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=12, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_x_discrete(labels=c("AM"="America", "EA"="Eurasia"))

svl.states <- ggplot(dt74, aes(x=state, y=log_svl)) + 
  stat_boxplot(geom = "errorbar", width = 0.25, lwd=0.6) +
  geom_boxplot(color="black", fill=pal, alpha=1, lwd=0.6, varwidth = T)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) 

svl.ecotype <- ggplot(dt74, aes(x=ecotype, y=log_svl)) + 
  geom_boxplot(color=c("dodgerblue4", "cornflowerblue", "chocolate4"), 
               fill=c("dodgerblue4", "cornflowerblue", "chocolate4"), alpha=0.6)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=11.5, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=25), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  ggtitle("Global")

svl.biome <-  ggplot(dt74, aes(x=biome, y=log_svl)) + 
  geom_boxplot(color=pal2, 
               fill=pal2, alpha=0.6)+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=11.5, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=25), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("Boreal forest"="Boreal","Desert"="Desert",
                            "Grasslands"="Grasslands", "Temperate Forest"="Temperate", 
                            "Tropical forest"="Tropical"))


library(cowplot)
plot_grid(
  svl.rad, svl.states, svl.ecotype, svl.biome,
  labels = c('A', 'B', 'C', 'D'),
  align="hv"
)

plot_grid(svl.rad, svl.states,
          labels=c("A", "B"),
          align = "hv",
          ncol=2,
          nrow=1)



# 3.2: Eurasia
#==============
svl.state.ea <- ggplot(dtEA, aes(x=state, y=log_svl)) + 
  geom_boxplot(color=c("darkorchid1", "darkorchid4"), fill=c("darkorchid1", "darkorchid4"), alpha=0.6)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=12, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_x_discrete(labels=c("AS"="Asia", "EU"="Europa"))+
  ggtitle("Eurasia")

svl.eco.ea <- ggplot(dtEA, aes(x=ecotype, y=log_svl)) + 
  geom_boxplot(color=c("dodgerblue4", "cornflowerblue", "chocolate4"), 
               fill=c("dodgerblue4", "cornflowerblue", "chocolate4"), alpha=0.6)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=11.5, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=25), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  ggtitle("Eurasia")+
  ylim(c(4,5))

svl.biome.ea <-  ggplot(dtEA, aes(x=biome, y=log_svl)) + 
  geom_boxplot(color=pal2[-2], 
               fill=pal2[-2], alpha=0.6)+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=11.5, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=25), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("Boreal forest"="Boreal",
                            "Grasslands"="Grasslands", "Temperate Forest"="Temperate", 
                            "Tropical forest"="Tropical"))+
  ylim(c(4,5))


library(cowplot)
plot_grid(
  svl.state.ea, svl.eco.ea, svl.biome.ea,
  labels = c('A', 'B', 'C'),
  align="hv"
)


# 3.3: America
#==============
svl.state.am <- ggplot(dtAM, aes(x=state, y=log_svl)) + 
  geom_boxplot(color=c("gold4", "gold3", "gold"), fill=c("gold4", "gold3", "gold"), alpha=0.6)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_x_discrete(labels=c("NAM"="North America", "MEX"="Mexico", "NEO"="Neotropics"))+
  ggtitle("America")

svl.eco.am <- ggplot(dtAM, aes(x=ecotype, y=log_svl)) + 
  geom_boxplot(color=c("dodgerblue4", "cornflowerblue", "chocolate4"), 
               fill=c("dodgerblue4", "cornflowerblue", "chocolate4"), alpha=0.6)+
  theme_classic() +
  ylab("Log(svl)")+
  xlab("") +
  theme(axis.text.x = element_text(size=11.5, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=25), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) + 
  ggtitle("America")

svl.biome.am <-  ggplot(dtAM, aes(x=biome, y=log_svl)) + 
  geom_boxplot(color=pal2, 
               fill=pal2, alpha=0.6)+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=11.5, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=25), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("Boreal forest"="Boreal",
                            "Grasslands"="Grasslands", "Temperate Forest"="Temperate", 
                            "Tropical forest"="Tropical"))


library(cowplot)
plot_grid(
  svl.state.ea, svl.eco.ea, svl.biome.ea,
  labels = c('A', 'B', 'C'),
  align="hv"
)


###########
# BIG PLOT
###########
# SVL differences by states, biomes and ecotypes in 3 rows: General, Eurasia and America

# General
plot_grid(svl.ecotype, svl.biome,
  # labels = c('A', 'B', 'C'),
  align="h", nrow=1, ncol=2
)

# Eurasia
plot_grid(svl.eco.ea, svl.biome.ea,
  # labels = c('A', 'B', 'C'),
  align="h", nrow=1, ncol=2
)

# America
plot_grid(svl.eco.am, svl.biome.am,
  # labels = c('A', 'B', 'C'),
  align="h", nrow=1, ncol=2
)



# 4. GGPLOT: Morphological variables
#=====================================
dt55 <- read.csv("DATA/dt_morpho_55sp.csv")
# Extract 6 species from EA radiation that are distributed in aM
dt55 <- dt55[dt55$sp!="Rana_tlaloci",]
dt55 <- dt55[dt55$sp!="Rana_luteiventris",]
dt55 <- dt55[dt55$sp!="Rana_boylii",]
dt55 <- dt55[dt55$sp!="Rana_sierrae",]
dt55 <- dt55[dt55$sp!="Rana_muscosa",]
dt55 <- dt55[dt55$sp!="Rana_aurora",]
dt55 <- dt55[dt55$sp!="Rana_cascadae",]
rownames(dt55) <- dt55$sp

library(reshape2)
palet <- c("gold3", "darkorchid")
morph <- melt(dt55[,c(4,14:17)], id = "phylo_reg") 
ggplot(morph, aes(x = variable, y = value, colour=phylo_reg, fill=phylo_reg)) +  # ggplot function
  geom_boxplot(alpha=0.5) +
  scale_colour_manual(values =palet, aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=13, face="bold"),
        axis.text.y = element_text(size=13))


library(reshape2)
morph.eco <- melt(dt55[,c(6,14:17)], id = "ecotype") 
ggplot(morph.eco, aes(x = variable, y = value, fill=ecotype, colour=ecotype)) +  # ggplot function
  geom_boxplot(alpha=0.7) +
  scale_colour_manual(values=c("dodgerblue4", "cornflowerblue", "chocolate4"),
                      aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10))


morph.state <- melt(dt55[,c(5,14:17)], id = "state") 
ggplot(morph.state, aes(x = variable, y = value, fill=state, colour=state)) +  # ggplot function
  geom_boxplot(alpha=0.6, position = position_dodge2(), width = 0.65) +
  scale_fill_manual(values=pal, aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10))

# Import biomes
biome <- read.csv("DATA/biome_68sp.csv")
rownames(biome) <- as.factor(biome$sp_lin)
biome <- biome[dt55$sp,]
biome2 <- na.omit(biome)
dt55$sp == biome2$sp_lin
dt55$biome <- biome2$new_biome

morph.biome <- melt(dt55[,c(14:17, 19)], id = "biome") 
ggplot(morph.biome, aes(x = variable, y = value, fill=biome, colour=biome)) +  # ggplot function
  geom_boxplot(alpha=0.6, position = position_dodge2(), width = 0.65) +
  scale_fill_manual(values=pal2, aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10)) +
  ggtitle("Global")

# Eurasia
morph.ea <- dt55[dt55$phylo_reg=="EA",]
morph.biome.ea <- melt(morph.ea[,c(14:17, 19)], id = "biome") 
ggplot(morph.biome.ea, aes(x = variable, y = value, fill=biome, colour=biome)) +  # ggplot function
  geom_boxplot(alpha=0.6, position = position_dodge2(), width = 0.65) +
  scale_fill_manual(values=pal2[-2], aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10)) +
  ggtitle("Eurasia")

# America 
morph.am <- dt55[dt55$phylo_reg=="AM",]
morph.biome.am <- melt(morph.am[,c(14:17, 19)], id = "biome") 
ggplot(morph.biome.am, aes(x = variable, y = value, fill=biome, colour=biome)) +  # ggplot function
  geom_boxplot(alpha=0.6, position = position_dodge2(), width = 0.6) +
  scale_fill_manual(values=pal2[-1], aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10)) +
  ggtitle("America")




# ADD subgenera 
subgen <- read.csv("DATA/subgen.csv", sep=";")
rownames(subgen) <- subgen$sp
subgen <- subgen[subgen$sp!="Rana_weiningensis",]
subgen <- subgen[subgen$sp!="Rana_shuchinae",]
subgen <- subgen[subgen$Subgenera!="",]
dt55$sp == subgen$sp

dt <- dt55[dt55$sp!="Rana_weiningensis",]
dt <- dt[dt$sp!="Rana_shuchinae",]

subgen <- na.omit(subgen[dt$sp,])
dt$sp == subgen$sp

dt$subgen <- subgen$Subgenera
pal3 <- c("#2D3142", "#4F5D75", "#BFC0C0", "#9E0031", "#EF8354")

morph.subgen <- melt(dt[,c(14:17,19)], id = "subgen") 
ggplot(morph.subgen, aes(x = variable, y = value, fill=subgen, colour=subgen)) +  # ggplot function
  geom_boxplot(alpha=0.7, position = position_dodge2()) +
  scale_fill_manual(values=pal3, aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10))

# SVL 
svl.subgen <- dt[,c(13,19),]
ggplot(svl.subgen, aes(x = subgen, y = log_svl, fill=subgen, colour=subgen)) +  # ggplot function
  geom_boxplot(alpha=0.7, position = position_dodge2()) +
  scale_fill_manual(values=pal3, aesthetics = c("colour", "fill"))+
  theme_classic() +
  ylab("Log(SVL)")+
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.title.y=element_text(size=13, face="bold"),
        axis.text.y = element_text(size=10))
  