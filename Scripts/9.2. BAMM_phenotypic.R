#############################
# BAMM: PHENOTYPIC EVOLUTION 
#############################

rm(list=ls())

library(phytools)
library(ape)
library(BAMMtools)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/BAMM/bamm-2.5.0-Windows")


#================
# 1. SVL (74 sp)
#================
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_68.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.5 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both ESS > 200

# summarize the posterior distribution of the number of shifts 
edata_trait <- getEventData("phyloRana_68sp.tre", eventdata = "event_data_trait_68.txt", 
                            burnin=0.5, type = "trait")

postfile <- "mcmc_out_trait_68.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)
layout(1,1)
plotPrior(postfile, expectedNumberOfShifts=1) # plot prior


# PLOT 1: phylorate plot
plot.bammdata(edata_trait, lwd=3, method="polar", pal="temperature")

# PLOT 2: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_trait, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 2218 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 3: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_trait, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 4: phenotypic evolutionary rate through time
tr68 <- read.tree("phyloRana_68sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr68))
plotRateThroughTime(edata_trait, intervalCol="black", 
                    avgCol="black",start.time=st, ylim=c(0,0.1), cex.axis=1)
text(x=30, y= 0.08, label="Rana (68 sp)", font=4, cex=1.5, pos=4)

# Separate radiations
#====================
allrates <- getCladeRates(edata_trait)
mean(allrates$beta) # tasa fenotípica = 0.008124642

plot(tr68, cex=0.7)
nodelabels(cex=0.5)

# Extraer MRCA americanas
american <- getCladeRates(edata_trait, node=100)
mean(american$beta)  # Mayor tasa fenotípica que todo el conjunto (0.0116)
mrca <- 100

# # Plot americanas
# plotRateThroughTime(edata_trait, intervalCol="darkgoldenrod", 
#                     avgCol="darkgoldenrod", xlim = c(st,0), ylim=c(0,0.1), cex.axis=1,
#                     node = mrca, nodetype="include")
# text(x=30, y= 0.08, label="SVL American (37 sp)", font=4, cex=1, pos=4)
# 
# # Plot Euroasiaticas
# plotRateThroughTime(edata_trait, intervalCol="darkorchid", 
#                     avgCol="darkorchid", xlim = c(st,0), ylim=c(0,0.1), cex.axis=1,
#                     node = mrca, nodetype="exclude")
# text(x=30, y= 0.08, label="SVL Eurasia (31 sp)", font=4, cex=1, pos=4)


# PLOT TOGETHER
layout(matrix(c(1,2,3), nrow=1, ncol=3))
plotRateThroughTime(edata_trait, intervalCol="black", 
                    avgCol="black",start.time=st, ylim=c(0,0.05), cex.axis=1)

plotRateThroughTime(edata_trait, intervalCol="darkorchid", 
                    avgCol="darkorchid", xlim = c(st,0), ylim=c(0,0.05), cex.axis=1,
                    node = mrca, nodetype="exclude")

plotRateThroughTime(edata_trait, intervalCol="gold", 
                    avgCol="gold", xlim = c(st,0), ylim=c(0,0.05), cex.axis=1,
                    node = mrca, nodetype="include")


# layout(1,1)

#==============
# 2. HW (51 sp)
#==============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_hw.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# Both ESS > 200

# summarize the posterior distribution of the number of shifts 
tr.hw <-read.tree("tree_hw.tre")
edata_hw<- getEventData(tr.hw, eventdata = "event_data_trait_hw.txt", 
                        burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_hw, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 35 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_hw, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr.hw))
plotRateThroughTime(edata_hw, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="HW Rana (51 sp)", font=4, cex=1, pos=4)


# Separate radiations
#====================
allrates <- getCladeRates(edata_hw)
mean(allrates$beta) # tasa fenotípica = 0.003064517

plot(tr.hw, cex=0.7)
nodelabels(cex=0.5)

# Extraer MRCA americanas
american <- getCladeRates(edata_hw, node=76)
mean(american$beta)  # tasa fenotípica = 0.00363984
mrca <- 76

# PLOT ALL together
layout(matrix(c(1,2,3), nrow=1, ncol=3))
plotRateThroughTime(edata_hw, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.01), cex.axis=1,
                    axis.labels = F)
plotRateThroughTime(edata_hw, intervalCol="darkorchid", 
                    avgCol="darkorchid", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="exclude", axis.labels = F)
plotRateThroughTime(edata_hw, intervalCol="gold",
                    avgCol="gold", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="include", axis.labels = F)


#==============
# 3. HL (45 sp)
#===============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_hl.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# Both ESS > 200

# summarize the posterior distribution of the number of shifts 
tr.hl<-read.tree("tree_hl.tre")
edata_hl<- getEventData(tr.hl, eventdata = "event_data_trait_hl.txt", 
                        burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_hl, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 7 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_hl, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr.hl))
plotRateThroughTime(edata_hl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="HL global (45 sp)", font=4, cex=1, pos=4)


# Separate radiations
#====================
allrates <- getCladeRates(edata_hl)
mean(allrates$beta) 

plot(tr.hl, cex=0.7)
nodelabels(cex=0.5)

# Extraer MRCA americanas
american <- getCladeRates(edata_hl, node=48)
mean(american$beta)  # tasa fenotípica = 0.00363984
mrca <- 48


# PLOT ALL together
layout(matrix(c(1,2,3), nrow=1, ncol=3))
plotRateThroughTime(edata_hl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.01), cex.axis=1,
                    axis.labels = F)
plotRateThroughTime(edata_hl, intervalCol="darkorchid", 
                    avgCol="darkorchid", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="exclude", axis.labels = F)
plotRateThroughTime(edata_hl, intervalCol="gold",
                    avgCol="gold", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="include", axis.labels = F)



#==============
# 4. FL (43 sp)
#==============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_fl.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# Both ESS > 200

# summarize the posterior distribution of the number of shifts 
tr.fl<-read.tree("tree_fl.tre")
edata_fl<- getEventData(tr.fl, eventdata = "event_data_trait_fl.txt", 
                        burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_fl, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 1 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_hl, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr.fl))
plotRateThroughTime(edata_fl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="FL global (43 sp)", font=4, cex=1, pos=4)


# Separate radiations
#====================
allrates <- getCladeRates(edata_fl)
mean(allrates$beta) 

plot(tr.fl, cex=0.7)
nodelabels(cex=0.5)

# Extraer MRCA americanas
american <- getCladeRates(edata_fl, node=46)
mean(american$beta)  # tasa fenotípica = 0.00363984
mrca <- 46


# PLOT ALL together
layout(matrix(c(1,2,3), nrow=1, ncol=3))
plotRateThroughTime(edata_fl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.01), cex.axis=1,
                    axis.labels = F)
plotRateThroughTime(edata_fl, intervalCol="darkorchid", 
                    avgCol="darkorchid", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="exclude", axis.labels = F)
plotRateThroughTime(edata_fl, intervalCol="gold",
                    avgCol="gold", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="include", axis.labels = F)



#==============
# 5. TL (46 sp)
#==============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_tl.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# Both ESS > 200

# summarize the posterior distribution of the number of shifts 
tr.tl<-read.tree("tree_tl.tre")
edata_tl<- getEventData(tr.tl, eventdata = "event_data_trait_tl.txt", 
                        burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_tl, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 1 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_tl, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr.tl))
plotRateThroughTime(edata_tl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="TL global (46 sp)", font=4, cex=1, pos=4)


# Separate radiations
#====================
allrates <- getCladeRates(edata_tl)
mean(allrates$beta) 

plot(tr.tl, cex=0.7)
nodelabels(cex=0.5)

# Extraer MRCA americanas
american <- getCladeRates(edata_tl, node=49)
mean(american$beta)  # tasa fenotípica = 0.00363984
mrca <- 49


# PLOT ALL together
layout(matrix(c(1,2,3), nrow=1, ncol=3))
plotRateThroughTime(edata_tl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.01), cex.axis=1,
                    axis.labels = F)
plotRateThroughTime(edata_tl, intervalCol="darkorchid", 
                    avgCol="darkorchid", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="exclude", axis.labels = F)
plotRateThroughTime(edata_tl, intervalCol="gold",
                    avgCol="gold", xlim = c(st,0), ylim=c(0,0.01), cex.axis=1,
                    node = mrca, nodetype="include", axis.labels = F)




##################
# PLOT ridgelines
####################
library("ggridges")
library("ggplot2")
div <- getTipRates(edata_trait)
diversity_data<- div$beta.avg
diversity_data <- as.data.frame(diversity_data)
diversity_data$group <- rownames(diversity_data)
colnames(diversity_data) <-  c("rate", "sp")

#Import traits 
dt68 <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/dt_68sp.csv")
diversity_data$sp == dt68$sp
diversity_data <- diversity_data[dt68$sp,]
diversity_data$radiation <- dt68$phylo_reg
diversity_data$state <- dt68$state
diversity_data$ecotype <- as.factor(dt68$ecotype)
diversity_data$biome <- as.factor(dt68$biome)

par(mar=c(6,6,6,6))
rad <- ggplot(diversity_data, aes(y = radiation, x = rate, fill=radiation)) +
  geom_density_ridges(scale=3, alpha=0.7)  +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_fill_cyclical(values = c("gold", "darkorchid")) +
  xlab("Rates of body size evolution") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

stat <- ggplot(diversity_data, aes(y = state, x = rate, fill=state)) +
  geom_density_ridges(scale=1, alpha=0.7)  +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_fill_cyclical(values = c("darkorchid1", "darkorchid4", "gold4", "gold3", "gold")) +
  xlab("Rates of body size evolution") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

eco <- ggplot(diversity_data, aes(y = ecotype, x = rate, fill=ecotype)) +
  geom_density_ridges(scale=0.9, alpha=0.7)  +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_fill_cyclical(values = c("dodgerblue4", "cornflowerblue", "chocolate4")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

biome <- ggplot(diversity_data, aes(y = biome, x = rate, fill=biome)) +
  geom_density_ridges(scale=0.9, alpha=0.7)  +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_fill_cyclical(values = c("cadetblue3", "antiquewhite3", "chocolate", "darkolivegreen4", "chartreuse3")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

library(cowplot)
plot_grid(
  rad, stat, eco, biome,
  labels = c('A', 'B', 'C', 'D'),
  align="hv"
)

subgen <- read.csv("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/DATA/subgen.csv", sep=";")
rownames(subgen) <- subgen$sp
subgen$sp == diversity_data$sp
subgen <- subgen[diversity_data$sp,]

diversity_data$subgen <- subgen$Subgenera
div2 <- diversity_data[diversity_data$sp!="Rana_weiningensis",]
div2 <- div2[div2$sp!="Rana_shuchinae",]
div2 <- div2[div2$subgen!="",]

pal <- c("#2D3142", "#4F5D75", "#BFC0C0", "#9E0031", "#EF8354")

subg <- ggplot(div2, aes(y = subgen, x = rate, fill=subgen)) +
  geom_density_ridges(scale=0.9, alpha=0.7)  +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_fill_cyclical(values = pal) +
  theme_bw() +
  ylab("") +
  xlab("Rates of body size evolution") +
  theme(axis.text.x = element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y = element_text(size=15), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
  
library(cowplot)
plot_grid(
  subg, stat, eco, biome,
  labels = c('A', 'B', 'C', 'D'),
  align="hv"
)
  

best <- getBestShiftConfiguration(edata_trait, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
library(RColorBrewer)
plot(best, lwd=2, color=c("red", "blue"))
addBAMMshifts(best, cex=1.5, col="white", bg="red")

par(mar=c(4,4,4,4))
plot.bammdata(edata_trait, lwd=3, method="phylogram", pal="RdBu", legend=T)
addBAMMshifts(best, cex=1.5, col="white", bg="red")

div2$rana <- as.factor("Rana")
ggplot(div2, aes(x = rate, y=rana, fill=stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_gradientn(colours = c("#053061", "#D1E5F0", "white","#F4A582", "#67001F"),
                      values = scales::rescale(c(0.0042, 0.02, 0.035, 0.051))) +
  theme_classic() + 
  ylab("Density") + 
  xlab("Distribution of body size evolutionary rate") +
  theme(axis.text.x = element_text(size=15), axis.title.x=element_text(size=20),
        axis.title.y = element_text(size=20))



