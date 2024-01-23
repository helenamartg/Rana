## NOTE: THIS IS AN UPDATED VERSION OF THE FUNCTION THAT IS COMPATIBLE WITH CHANGES MADE IN GEIGER 1.99 ##

## Written by Graham Slater 10/ 10 / 2013 ##


## This function computes the average subclade disparity through time and plots it, a la Harmon et al. 2003, but adds the 95% range for the simulations shaded by "color" and a P value returned a la Slater et al. 2010. I also modified the plot so that the median, rather than the mean of the simulations is plotted as the dashed line.


dttFullCIs<-function (phy, data, disp = "avg.sq", nsims = 1000, mdi.range = c(0, 1),ciRange=c(0.025,0.975), color = "gray81", linecol="black") 
{
    td <- treedata(phy, data, sort = T);
    dtt.data <- dtt(td$phy, td$data, index = disp, plot = F)
    # dtt.data2 <- as.data.frame(dtt.data$times)
    # dtt.data2$branches <- rownames(dtt.data2)
    # branch.times<- as.data.frame(c(0, branching.times(td$phy)))
    # names(branch.times) <- "times"
    # dtt.data2$new_times <- branch.times$times
    # dtt.data3 <- as.matrix(dtt.data2)
    # dtt.data4 <- dtt.data3[,3]
    # dtt.data4 <- as.factor(dtt.data4)
    # dtt.data4 <- as.data.frame(dtt.data4)
    # # 
    #  
    # dim(branch.times) == dim(dtt.data$times)
    # new_times <- dtt.data$times
    # new_times$nodes <- rownames(new_times)
    # branch.times$nodes <- rownames(branch.times)
    # branch.times$nodes==new_times$nodes
    # new_times <- branch.times[new_times$nodes,]
    # new_times <- new_times[-1,] 
    # new_times <- rbind(0, new_times)
    # dtt.data$times <- dtt.data4
    
   plot(dtt.data$times, dtt.data$dtt, type = "l",col="white",lwd = 2, xlab = "Relative time", ylab = "Disparity", ylim=c(0,2))
    s<-vcv(td$phy, td$dat)
    sims <- sim.char(td$phy, s, nsims)
    dtt.sims <- sapply(apply(sims,3, dtt, phy=td$phy, index=disp, plot=F), "[[", "dtt");
    
    
    quantiles<-numeric()
for (i in 1:length(dtt.data$times)){
	
	range<-quantile(dtt.sims[i,],ciRange)
	quantiles<-rbind(quantiles,range)
	}
	
median.sims <- apply(dtt.sims, 1, median)
	 
x<-rep(dtt.data$times,2)
y<-(quantiles[,1])
y<-append(y,quantiles[,2])
polygon(x,y,border=color,col = color)

lines(dtt.data$times, median.sims, lty = 2,lwd=1.5)
lines(dtt.data$times, dtt.data$dtt,lty= 1,lwd=2.5, col=linecol)

    MDI <- geiger:::.area.between.curves(dtt.data$times, apply(dtt.sims, 1, median),dtt.data$dtt, mdi.range)
        
 as.data.frame(dtt.sims)->simframe

mdiNull<-numeric(nsims)
for (i in 1:nsims){
	MDIsim <- geiger:::.area.between.curves(dtt.data$times, simframe[,i],dtt.data$dtt, mdi.range)
	mdiNull[i]<-MDIsim
	}

if(MDI<0){
which(mdiNull>0)->pos
(length(pos)/nsims)->prop

}
else {which(mdiNull<0)->pos
(length(pos)/nsims)->prop
}
    return(list(dtt.data = dtt.data, dtt.sims = dtt.sims, times = dtt.data$times, MDI.nullDistribution = mdiNull,
        MDI = MDI,Pvalue=prop))
}

