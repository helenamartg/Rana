#############################
# BAMM: PHENOTYPIC EVOLUTION 
#############################

rm(list=ls())

library(phytools)
library(ape)
library(BAMMtools)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/BAMM/bamm-2.5.0-Windows")


##########
# 1. SVL
##########

# 1.1 All species (74sp)
#=======================
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_200.txt", header=T)
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
edata_trait <- getEventData("phyloRana_74sp.tre", eventdata = "event_data_trait_200.txt", 
                            burnin=0.5, type = "trait")

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
tr74 <- read.tree("phyloRana_74sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr74))
plotRateThroughTime(edata_trait, intervalCol="black", 
                    avgCol="black",start.time=st, ylim=c(0,0.1), cex.axis=1)
text(x=30, y= 0.08, label="Rana (74 sp)", font=4, cex=1.5, pos=4)

# Separate radiations
#====================
allrates <- getCladeRates(edata_trait)
mean(allrates$beta) # tasa fenotípica = 0.0088

plot(tr74, cex=0.7)
nodelabels()

# Extraer MRCA americanas
american <- getCladeRates(edata_trait, node=112)
mean(american$beta)  # Mayor tasa fenotípica que todo el conjunto (0.0116)
mrca <- 112

# Plot americanas
plotRateThroughTime(edata_trait, intervalCol="darkgoldenrod", 
                    avgCol="darkgoldenrod", xlim = c(st,0), ylim=c(0,0.1), cex.axis=1,
                    node = mrca, nodetype="include")
text(x=30, y= 0.08, label="American (37 sp)", font=4, cex=1.5, pos=4)

# Plot Euroasiaticas
plotRateThroughTime(edata_trait, intervalCol="darkorchid", 
                    avgCol="darkorchid", xlim = (st,0), ylim=c(0,0.1), cex.axis=1,
                    node = mrca, nodetype="exclude")
text(x=30, y= 0.08, label="Eurasia (31 + 6sp)", font=4, cex=1.5, pos=4)


# 1.2  EURASIA 
#==============
#Assessing convergence
mcmcout <- read.csv("mcmc_out_trait_EA_200.txt", header=T)
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
trEA <- read.tree("phyloEA_31sp.tre")
edata_trait_EA <- getEventData(trEA, eventdata = "event_data_trait_EA_200.txt", 
                               burnin=0.5, type = "trait")
shift_probs <- summary(edata_trait_EA)

# PLOT 2: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_trait_EA, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 3 distinct shift configurations
summary(cset)
plot.credibleshiftset(cset, lwd=2.5)

# PLOT 3: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_trait_EA, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 4: phenotypic evolutionary  rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(trEA))
plotRateThroughTime(edata_trait_EA, intervalCol="orchid", 
                    avgCol="orchid", start.time=st, ylim=c(0,0.1), cex.axis=1)
text(x=30, y= 0.08, label="Eurasia (31 sp)", font=4, cex=1.5, pos=4)



# 1.3 AMERICA 
#=============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_trait_AM_200.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both ESS > 200

# summarize the posterior distribution of the number of shifts 
trAM <-read.tree("phyloAM_37sp.tre")
edata_trait_AM <- getEventData(trAM, eventdata = "event_data_trait_AM_200.txt", 
                               burnin=0.5, type = "trait")
shift_probs <- summary(edata_trait_AM)

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_trait_AM, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 290 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_trait_AM, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(trAM))
plotRateThroughTime(edata_trait_AM, intervalCol="gold", 
                    avgCol="gold", start.time=st, ylim=c(0,0.1), cex.axis=1)
text(x=30, y= 0.08, label="America (37 sp)", font=4, cex=1.5, pos=4)



################
# 2. HW 
###############

# 2.1 All species (55sp)
#================
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_hw_200.txt", header=T)
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
tr55 <-read.tree("phyloRana_55sp.tre")
edata_hw<- getEventData(tr55, eventdata = "event_data_hw_200.txt", 
                        burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_hw, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 25 distinct shift configurations
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
st <- max(branching.times(tr55))
plotRateThroughTime(edata_hw, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="HW global (55 sp)", font=4, cex=1.5, pos=4)



# 2.2 EURASIA
#============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_hw_EA_200.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both ESS > 200

# summarize the posterior distribution of the number of shifts 
phyEA <-read.tree("phyEA_24sp.tre")
edata_hw_EA <- getEventData(phyEA, eventdata = "event_data_hw_EA_200.txt", 
                            burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_hw_EA, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 3 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_hw_EA, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(phyEA))
plotRateThroughTime(edata_hw_EA, intervalCol="darkorchid", 
                    avgCol="darkorchid", start.time=st, ylim=c(0,0.1), cex.axis=1)
text(x=30, y= 0.08, label="HW Eurasia \n (24 sp)", font=4, cex=1.5, pos=4)



# 2.3 AMERICA
#=============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_hw_AM_200.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 50% of burnin
burnstart <- floor(0.50 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both ESS > 200

# summarize the posterior distribution of the number of shifts 
phyAM <-read.tree("phyAM_27sp.tre")
edata_hw_AM <- getEventData(phyAM, eventdata = "event_data_hw_AM_200.txt", 
                            burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_hw_AM, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 8 distinct shift configurations
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_hw_AM, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: phenotypic evolutionary rate through time
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(phyAM))
plotRateThroughTime(edata_hw_AM, intervalCol="gold", 
                    avgCol="gold", start.time=st, ylim=c(0,0.1), cex.axis=1)
text(x=30, y= 0.08, label="HW America (27 sp)", font=4, cex=1.5, pos=4)




###########
# 3. HL
###########

# 3.1 All species (48sp)
#===============
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_hl_200.txt", header=T)
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
tr48 <-read.tree("tr48_hl.tre")
edata_hl<- getEventData(tr48, eventdata = "event_data_hl_200.txt", 
                        burnin=0.5, type = "trait")

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_hl, expectedNumberOfShifts=1, threshold=3)
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
st <- max(branching.times(tr48))
plotRateThroughTime(edata_hl, intervalCol="black", 
                    avgCol="black", start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="HL global (48 sp)", font=4, cex=1.5, pos=4)


# 3.2 EURASIA
#==============


# 3.3 AMERICA
#===============
plot(tr48, cex=0.5)
nodelabels()

rates_hl <- getCladeRates(edata_hl)

mean(rates_hl$beta)
quantile(rates_hl$beta, c(0.05, 0.95))

# Extraer tasas de sp americanas
american <- getCladeRates(edata_hl, node=76)
mean(american$beta)

# Seleccionar dos species americanas
species1 <- "Rana_yavapaiensis"
species2 <- "Rana_okaloosae"

tipnode1 <- which(tr48$tip.label == species1)
tipnode2 <- which(tr48$tip.label == species2)

# Extraer mrca de las 2 sp
mrca <- getMRCA(tr48, tip = c(tipnode1, tipnode2))

# Plotear rate through time de las americanas
plotRateThroughTime(edata_hl, node = mrca, nodetype="include", ylim=c(0,0.7))
text(x=20, y=0.5, "HL American", font=4, cex=1.5, pos=4)

# 3.4 EUROPEAS
#==============
plotRateThroughTime(edata_hl, node = mrca, nodetype="exclude", ylim=c(0,0.7))
text(x=20, y=0.5, "HL European \n (including 6sp)", font=4, cex=1.5, pos=4)





