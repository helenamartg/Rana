#######################################
# BAMM: DIVERSIFICATION
#######################################

rm(list=ls())

library(phytools)
library(ape)
library(BAMMtools)

setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/Rana_genus_data/GITHUB/BAMM/bamm-2.5.0-Windows")


#======================
# 1. All species (74sp)
#=======================
# Assessing convergence (RUN 2 with 100 million generations)
mcmcout <- read.csv("mcmc_out2.txt", header=T)
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
edata <- getEventData("phyloRana_74sp.tre", eventdata = "event_data2.txt", burnin=0.25)
shift_probs <- summary(edata)

# prior distribution in BAMM --> bayes factor
postfile <- "mcmc_out2.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)

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
tr74 <- read.tree("phyloRana_74sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(tr74))
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st, 
                    ylim=c(0,0.5), cex.axis=1, ratetype = "netdiv")
text(x=30, y= 0.4, label="Rana (74 sp)", font=4, cex=2.0, pos=4)


# Separate by radiation
#=======================
# AMERICA
allrates <- getCladeRates(edata)
mean(allrates$lambda) # 0.07037

plot(tr74, cex=0.7)
nodelabels()  # Node 112 = american sp

# Extraer MRCA americanas
american <- getCladeRates(edata, node=112)
mean(american$lambda) # 0.0691

mrca <- 112

# Plot americans
plotRateThroughTime(edata, intervalCol="darkgoldenrod", avgCol="darkgoldenrod", 
                    node = mrca, nodetype="include", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
text(x=35, y= 0.4, label="American (37 sp)", font=4, cex=2.0, pos=4)

# Plot europeans
plotRateThroughTime(edata, intervalCol="darkorchid", avgCol="darkorchid",
                    node = mrca, nodetype="exclude", ylim=c(0,0.5), 
                    xlim=c(st,0), ratetype = "netdiv")
text(x=25, y=0.4, "European \n (including 6sp)", font=4, cex=1.5, pos=4)




#===================================
# 2. Diversification EURASIA (31 sp)
#===================================
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_div_EA.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 10% of burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both ESS > 200

# summarize the posterior distribution of the number of shifts 
edata_div_EA <- getEventData("phyloEA_31sp.tre", eventdata = "event_data_div_EA.txt", 
                            burnin=0.1)

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_div_EA, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 3 distinct shift configurations
summary(cset)
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_div_EA, expectedNumberOfShifts = 1,
                                  threshold = 5)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: diversification rate through time
trEA <- read.tree("phyloEA_31sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(trEA))
plotRateThroughTime(edata_div_EA, intervalCol="orchid", 
                    avgCol="orchid", ratetype = "netdiv",
                    start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="Eurasia (31 sp)", font=4, cex=1.5, pos=4)



#===================================
# 2. Diversification AMERICA (37 sp)
#===================================
#Assessing convergence 
mcmcout <- read.csv("mcmc_out_div_AM.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation, type="l")

#Discard 10% of burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation, type="l")

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both ESS > 200

# summarize the posterior distribution of the number of shifts 
edata_div_AM <- getEventData("phyloAM_37sp.tre", eventdata = "event_data_div_AM.txt", 
                             burnin=0.1)

# PLOT 1: set of shift configurations (that account for 95% of the prob) and their frequencies
cset <- credibleShiftSet(edata_div_AM, expectedNumberOfShifts=1, threshold=3)
cset$number.distinct  # 3 distinct shift configurations
summary(cset)
plot.credibleshiftset(cset, lwd=2.5, plotmax = 4)

# PLOT 2: Get the rate shift configuration with the maximum a posteriori probability, 
best <- getBestShiftConfiguration(edata_div_AM, expectedNumberOfShifts = 1,
                                  threshold = 5)
layout(1,1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)

# PLOT 3: diversification rate through time
trAM <- read.tree("phyloAM_37sp.tre")
par(mar=c(6,6,6,6))
layout(1,1)
st <- max(branching.times(trAM))
plotRateThroughTime(edata_div_AM, intervalCol="gold", 
                    avgCol="gold", ratetype = "netdiv",
                    start.time=st, ylim=c(0,0.5), cex.axis=1)
text(x=30, y= 0.4, label="America (37 sp)", font=4, cex=1.5, pos=4)








