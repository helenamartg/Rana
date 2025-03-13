#######################################
# BAMM: DIVERSIFICATION
#######################################

rm(list=ls())

library(phytools)
library(ape)
library(BAMMtools)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/BAMM/bamm-2.5.0-Windows")


#===========================================================
# 1. All species (84sp): Whole Yuan phylogeny without 6 sp
#===========================================================
# Assessing convergence (RUN 2 with 100 million generations)
mcmcout <- read.csv("mcmc_out_div_83.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

# Discard 25% of burnin
burnstart <- floor(0.25 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# ESS: Check the effective sample sizes of the log-likelihood 
#  and the number of shift events present in each sample
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# Both ESS > 200

# model posterior probabilities
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs) # 6 models

# summarize the posterior distribution of the number of shifts 
edata <- getEventData("phyloRana_83sp.tre", eventdata = "event_data_div_83.txt", burnin=0.25)
shift_probs <- summary(edata)

# prior distribution in BAMM --> bayes factor
postfile <- "mcmc_out_div_83.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)
layout(1,1)
plotPrior(postfile, expectedNumberOfShifts=1) # plot prior

marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs) # plot marginal odss for rate shift locations

# PLOT 1: phylorate plot
plot.bammdata(edata, lwd=3, method="polar", pal="temperature")

# PLOT 2: shift configurations and their frequencies
cset <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3)
plot.credibleshiftset(cset, lwd=2.5)
# set of distinct shift configurations that account for 95% of the probability of the data

# PLOT3: best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 4: diversification rate through time
tr83 <- read.tree("phyloRana_83sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr83))
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st, 
                    ylim=c(0,0.5), cex.axis=1, ratetype = "netdiv")
text(x=30, y= 0.4, label="Rana (83 sp)", font=4, cex=2.0, pos=4)


# Separate by radiation
#=======================
# AMERICA
allrates <- getCladeRates(edata)
mean(allrates$lambda) # 0.076824

par(mar=c(1,1,1,1))
plot(tr83, cex=0.7)
nodelabels()  # Node 86 = american sp

# Extraer MRCA americanas
american <- getCladeRates(edata, node=86)
mean(american$lambda) # 0.075347

mrca <- 86

# Plot americans
plotRateThroughTime(edata, intervalCol="darkgoldenrod", avgCol="darkgoldenrod", 
                    node = mrca, nodetype="include", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
text(x=35, y= 0.4, label="America", font=4, cex=2.0, pos=4)

# Plot europeans
plotRateThroughTime(edata, intervalCol="darkorchid", avgCol="darkorchid",
                    node = mrca, nodetype="exclude", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
text(x=25, y=0.4, "Eurasia", font=4, cex=2, pos=4)




#=====================================
# 2. Reduced (74sp) without 6 species
#=====================================
# Assessing convergence (RUN 2 with 100 million generations)
mcmcout <- read.csv("mcmc_out_div_68.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

# Discard 25% of burnin
burnstart <- floor(0.25 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# ESS: Check the effective sample sizes of the log-likelihood 
#  and the number of shift events present in each sample
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# Both ESS > 200

# model posterior probabilities
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs) # 5 models

# summarize the posterior distribution of the number of shifts 
edata <- getEventData("phyloRana_68sp.tre", eventdata = "event_data_div68.txt", burnin=0.25)
shift_probs <- summary(edata)

# prior distribution in BAMM --> bayes factor
postfile <- "mcmc_out_div_68.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)
layout(1,1)
plotPrior(postfile, expectedNumberOfShifts=1) # plot prior

marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs) # plot marginal odss for rate shift locations

# PLOT 1: phylorate plot
plot.bammdata(edata, lwd=3, method="polar", pal="temperature")

# PLOT 2: shift configurations and their frequencies
cset <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3)
plot.credibleshiftset(cset, lwd=2.5)
# set of distinct shift configurations that account for 95% of the probability of the data

# PLOT3: best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 1,
                                  threshold = 3)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 4: diversification rate through time
tr68 <- read.tree("phyloRana_68sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr68))
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st, 
                    ylim=c(0,0.5), cex.axis=1, ratetype = "netdiv")
text(x=30, y= 0.4, label="Rana (68 sp)", font=4, cex=2.0, pos=4)


# Separate by radiation
#=======================
# AMERICA
allrates <- getCladeRates(edata)
mean(allrates$lambda) # 0.07053139

par(mar=c(1,1,1,1))
plot(tr68, cex=0.7)
nodelabels()  # Node 100 = american sp

# Extraer MRCA americanas
american <- getCladeRates(edata, node=100)
mean(american$lambda) # 0.0691727

mrca <- 100

# Plot americans
plotRateThroughTime(edata, intervalCol="darkgoldenrod", avgCol="darkgoldenrod", 
                    node = mrca, nodetype="include", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
text(x=35, y= 0.4, label="America (37sp)", font=4, cex=2.0, pos=4)

# Plot europeans
plotRateThroughTime(edata, intervalCol="darkorchid", avgCol="darkorchid",
                    node = mrca, nodetype="exclude", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
text(x=25, y=0.4, "Eurasia (31 sp)", font=4, cex=2, pos=4)




# PLOT THREE
layout(matrix(c(1,2,3), nrow=1, ncol=3))

plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st, 
                    ylim=c(0,0.5), cex.axis=1, ratetype = "netdiv")

plotRateThroughTime(edata, intervalCol="darkgoldenrod", avgCol="darkgoldenrod", 
                    node = mrca, nodetype="include", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")

plotRateThroughTime(edata, intervalCol="darkorchid", avgCol="darkorchid",
                    node = mrca, nodetype="exclude", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
