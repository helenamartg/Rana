##############################
# CLeaning datasets
############################


rm(list=ls())
library(phytools)
library(ape)
library(castor)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB")

# PHYLOGENY
#-----------
phyloRana <- read.tree("DATA/phyloRana_74sp.tre")
phyloRana$tip.label

# There are 6 speciesin EA radiation that are actually in America
# after removing, there must remain 68 species
# Rana_luteiventris
# Rana_boylii
# Rana_sierrae
# Rana_muscosa
# Rana_aurora
# Rana_cascadae

phyloRana$tip.label
phyloRana_68sp <- drop.tip(phyloRana, c(3:8))
plot(phyloRana_68sp)

# SAVE tree
# write.tree(phyloRana_68sp, "phyloRana_68sp.tre")

phyloAM <- extract.clade(phyloRana_68sp, 100)
plot(phyloAM)
phyloEA <- drop.tip(phyloRana_68sp, c(32:68))
plot(phyloEA)

plot(phyloAM)
is.ultrametric(phyloAM)

plot(phyloEA)
is.ultrametric(phyloEA)

phyloRana$tip.label
morpho <- read.csv("DATA/dt_morpho_55sp.csv")
morpho <- morpho[morpho$sp!="Rana_tlaloci",]
morpho$sp
phylo55sp <- get_subtree_with_tips(phyloRana, only_tips = morpho$sp, force_keep_root = T)
plot(phylo55sp$subtree)

length(morpho$sp) == length(phylo55sp$subtree$tip.label)


# write.tree(phyloAM, "phyloAM_37sp.tre")
# write.tree(phyloEA, "phyloEA_31sp.tre")
# write.tree(phylo55sp$subtree, "phyloRana_55sp.tre")


# TRAIT DATASET
#---------------
dt74 <- read.csv("DATA/dt_74sp.csv", header=T)
rownames(dt74) <- dt74$sp
dt74 <- dt74[,-1]

dt68 <- dt74[dt74$sp!="Rana_luteiventris",]
dt68 <- dt68[dt68$sp!="Rana_boylii",]
dt68 <- dt68[dt68$sp!="Rana_sierrae",]
dt68 <- dt68[dt68$sp!="Rana_muscosa",]
dt68 <- dt68[dt68$sp!="Rana_aurora",]
dt68 <- dt68[dt68$sp!="Rana_cascadae",]

# BIOME
biome <- read.csv("DATA/raw_data/biome_all.csv")
biome <- biome[biome$sp_lin!="Rana_psilonota",]

rownames(biome) <- biome$sp_lin
biome$sp_lin == dt74$sp
biome <- biome[dt74$sp,]
dt74$biome <- biome$new_biome

# write.csv(dt74, "dt_74sp.csv")

biome <- biome[biome$sp!="Rana_luteiventris",]
biome <- biome[biome$sp!="Rana_boylii",]
biome <- biome[biome$sp!="Rana_sierrae",]
biome <- biome[biome$sp!="Rana_muscosa",]
biome <- biome[biome$sp!="Rana_aurora",]
biome <- biome[biome$sp!="Rana_cascadae",]

biome$sp_lin==dt68$sp
biome <- biome[dt68$sp,]
dt68$biome <- biome$new_biome

# SAVE
# write.csv(dt68, "dt_68sp.csv")



# CLIMATIC DATASET
clim_AM <- read.csv("DATA/env_vars_AM_sp.csv", sep=";")
rownames(clim_AM) <- clim_AM$X
clim_AM <- clim_AM[clim_AM$X!="Rana_psilonota",]

clim_EA <- read.csv("DATA/env_vars_EA_sp.csv", sep=";")
rownames(clim_EA) <- clim_EA$X
clim <- rbind(clim_AM, clim_EA)

clim_EA <- clim_EA[clim_EA$X!="Rana_luteiventris",]
clim_EA <- clim_EA[clim_EA$X!="Rana_boylii",]
clim_EA <- clim_EA[clim_EA$X!="Rana_sierrae",]
clim_EA <- clim_EA[clim_EA$X!="Rana_muscosa",]
clim_EA <- clim_EA[clim_EA$X!="Rana_aurora",]
clim_EA <- clim_EA[clim_EA$X!="Rana_cascadae",]



# SAVE
# write.csv(clim_AM, "clim_AM_37sp.csv")
# write.csv(clim_EA, "clim_EA_31sp.csv")
# write.csv(clim, "clim_74sp.csv")

