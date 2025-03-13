##########################
# Correlation tests
##########################

rm(list=ls())
setwd("C:/Users/helen/Desktop/MASTER_UB/TFM/4Publication/2nd_Revision/Reanalisis")

# 1. Import climatic variables
#===============================
clim.EA <- read.csv("DATA/climEA_26sp.csv")
clim.AM <- read.csv("DATA/climAM_33sp.csv")

clim <- rbind(clim.EA, clim.AM)
clim <- clim[,-1]
rownames(clim) <- clim$X


# 2. Scale climatic variables
#============================
# first remove SD
clim2 <- clim[,-grep("SD", colnames(clim))]
clim3 <- scale(clim2[,-1])
colnames(clim3)


# 3. Correlation matrix
#=======================
library(Hmisc)
results <- rcorr(as.matrix(clim3))

# Extract correlation coefficients
correlation_matrix <- results$r

# Extract p-values
p_value_matrix <- results$P

mat <- matrix(ncol=27, nrow=27)
rownames(mat) <- rownames(p_value_matrix)
colnames(mat) <- colnames(p_value_matrix)

# Substitute NAs od the diagonal
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (is.na(p_value_matrix[i,j])){
      p_value_matrix[i,j] <- 1
    } else {}
  }
}

# Only significant correlations
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (p_value_matrix[i,j] <= 0.05){
      mat[i,j] <- correlation_matrix[i,j]
      
    } else {
      mat[i,j] <- 0
    }
  }
}

# Only correlations with r > 0.76
for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    if (mat[i,j] <= 0.76){
      mat[i,j] <- 0
    }
  }
}

mat <- as.data.frame(mat)


# SAVE
# write.csv(correlation_matrix, "Results/Cor_matrix_climvars.csv")
# write.csv(p_value_matrix, "Results/Cor_pvalue_climvars.csv")
# write.csv(mat, "Results/Cor_onlysignif_climvars.csv")


# 1). Bergmann's rule # 
#=======================
# Solar radiation + Temperature + PET + ELEV
bergman <- clim3[,c(7:9,13:18, 22:24)]
results <- rcorr(as.matrix(bergman))
round(results$P,3) 
round(results$r, 3)

# Extract correlation coefficients
correlation_matrix <- results$r

# Extract p-values
p_value_matrix <- results$P

mat <- matrix(ncol=12, nrow=12)
rownames(mat) <- rownames(p_value_matrix)
colnames(mat) <- colnames(p_value_matrix)

# Substitute NAs od the diagonal
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (is.na(p_value_matrix[i,j])){
      p_value_matrix[i,j] <- 1
    } else {}
  }
}

# Only significant correlations
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (p_value_matrix[i,j] <= 0.05){
      mat[i,j] <- correlation_matrix[i,j]
      
    } else {
      mat[i,j] <- 0
    }
  }
}

# Only correlations with r > 0.76
for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    if (mat[i,j] <= 0.76){
      mat[i,j] <- 0
    }
  }
}

mat <- as.data.frame(mat)
plot(hclust(dist(mat)), main="Bergmann")


# selected variables: TMAXavg + TMINavg + SRADmax + SRADmin
bergman <- clim3[,c(14,15, 17, 18)]
results <- rcorr(as.matrix(bergman))
round(results$P,3) # All significant
round(results$r, 3) # Still correlated

plot(hclust(dist(results$r)), main="Bergmann")


# reselected variables: TMINavg + SRADavg + PETmax  + ELEVmax + ELEVmin
bergman <- clim3[,c(8,13,18,23,24)]
results <- rcorr(as.matrix(bergman))
round(results$P,3) # Significance
round(results$r, 3) # NOT CORRELATED



# 2). Water conservation # 
#========================
# VAPR + PREC + PÊT + ARID
water <- clim3[,c(1:3, 7:12, 19:21)]
results <- rcorr(as.matrix(water))
round(results$P,3) # Significance
round(results$r, 3) # Correlations

mat <- results$r
plot(hclust(dist(mat)), main="water conservation")

# Extract p-values
p_value_matrix <- results$P
correlation_matrix <- results$r

mat <- matrix(ncol=12, nrow=12)
rownames(mat) <- rownames(p_value_matrix)
colnames(mat) <- colnames(p_value_matrix)

# Substitute NAs of the diagonal
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (is.na(p_value_matrix[i,j])){
      p_value_matrix[i,j] <- 1
    } else {}
  }
}

# Only significant correlations
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (p_value_matrix[i,j] <= 0.05){
      mat[i,j] <- correlation_matrix[i,j]
      
    } else {
      mat[i,j] <- 0
    }
  }
}

# Only correlations with r > 0.76
for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    if (mat[i,j] <= 0.76){
      mat[i,j] <- 0
    }
  }
}

mat <- as.data.frame(mat)  # Significant correlated variables



# Selected variables: ARIDavg + ARIDmax + ARIDmin + VAPRmax + VAPRmin + 
#                     PETmin + PRECmax + PRECavg

water <- clim3[,c(1:3, 9:11, 20:21)]
results <- rcorr(as.matrix(water))
round(results$P,3) # Significance
round(results$r, 3) # Correlations


# Extract p-values
p_value_matrix <- results$P

mat <- matrix(ncol=8, nrow=8)
rownames(mat) <- rownames(p_value_matrix)
colnames(mat) <- colnames(p_value_matrix)

# Substitute NAs of the diagonal
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (is.na(p_value_matrix[i,j])){
      p_value_matrix[i,j] <- 1
    } else {}
  }
}


# Only significant correlations
for (i in 1:nrow(p_value_matrix)){
  for (j in 1:ncol(p_value_matrix)){
    
    if (p_value_matrix[i,j] <= 0.05){
      mat[i,j] <- correlation_matrix[i,j]
      
    } else {
      mat[i,j] <- 0
    }
  }
}

# Only correlations with r > 0.76
for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    if (mat[i,j] <= 0.76){
      mat[i,j] <- 0
    }
  }
}

mat <- as.data.frame(mat)  # NOT CORRELATED VARS 

plot(hclust(dist(mat)))




# 3). Primary productivity hypothesis # 
#=====================================
# NDVI
productivity <- clim3[,4:6]
results <- rcorr(as.matrix(productivity))
round(results$P,3) # All significant
round(results$r, 3) # Correlations

mat <- results$r
plot(hclust(dist(mat)), main="Productivity")

# Selected: NDVImax + NDVImin
productivity <- clim3[,5:6]
results <- rcorr(as.matrix(productivity))
round(results$P,3) # All significant
round(results$r, 3) # Correlations


# 4). Habitat availability #
#===========================
# ELEV
habitat <- clim3[,22:24]
results <- rcorr(as.matrix(habitat))
round(results$P,3) # Almost all significant
round(results$r, 3) # Correlations

plot(hclust(dist(results$r)))

# Selected: ELEVmax + ELEVmin
habitat <- clim3[,23:24]
results <- rcorr(as.matrix(habitat))
round(results$P,3) # Not significant
round(results$r, 3) # NOT CORRELATED




# 12 SELECTED VARIABLES non correlated
#========================================
# From all 27 variables we selected 12 non-correlated between them. Instead of testing by hypothesis
# Run correlation test between selected variables
# from 27 to 12! 

clim3 <- as.data.frame(clim3)
clim.red <- clim3[c("ARIDavg", "ARIDmin", "NDVImax", "PETmax", "PRECmax", "PRECmin",
                    "SRADavg", "TMINavg", "VAPRmax", "ELEVmax", "ELEVmin",
                    "WINDmin")]

results2 <- rcorr(as.matrix(clim.red))
correlation_matrix2 <- results2$r
p_value_matrix2 <- results2$P

correlation_matrix2 <- as.data.frame(correlation_matrix2)


mat <- matrix(ncol=12, nrow=12)
rownames(mat) <- rownames(p_value_matrix2)
colnames(mat) <- colnames(p_value_matrix2)

# Substitute NAs of the diagonal
for (i in 1:nrow(p_value_matrix2)){
  for (j in 1:ncol(p_value_matrix2)){
    
    if (is.na(p_value_matrix2[i,j])){
      p_value_matrix2[i,j] <- 1
    } else {}
  }
}


# Only significant correlations
for (i in 1:nrow(p_value_matrix2)){
  for (j in 1:ncol(p_value_matrix2)){
    
    if (p_value_matrix2[i,j] <= 0.05){
      mat[i,j] <- correlation_matrix2[i,j]
      
    } else {
      mat[i,j] <- 0
    }
  }
}

# Only correlations with r > 0.76
for (i in 1:nrow(mat)){
  for (j in 1:ncol(mat)){
    if (mat[i,j] <= 0.76){
      mat[i,j] <- 0
    }
  }
}

mat <- as.data.frame(mat)  # NOT CORRELATED VARS 

plot(hclust(dist(mat)))
# 12 no correlated vars



# Save non-correlated variables
p_value_matrix2 <- round(mat, 3)
# write.csv(p_value_matrix2, "Results/clim_vars_noncorrelated.csv")
