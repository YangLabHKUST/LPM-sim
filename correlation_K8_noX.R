##### Performance in characterizing the correlation among the traits #####
##### when annotations have no role (eight traits)                   #####
# Supplementary Figure S2

library(MASS)
library(pbivnorm)
library(mvtnorm)

K <- 8                 # No. of traits
M <- 100000            # No. of SNPs
beta0 <- -1            # intercept of the probit model
beta0 <- rep(beta0, K)
set.seed(1)

alpha <- c(0.2, 0.35, 0.5, 0.3, 0.45, 0.55, 0.25, 0.4) # parameter in the Beta distribution
R <- matrix(0, K, K)   # Correlation matrix for the traits
R[1, 2] <- 0.7
R[1, 3] <- 0.4
R[2, 3] <- 0.2
R[4, 5] <- 0.6
R[4, 6] <- 0.3
R[5, 6] <- 0.1
R[7, 8] <- 0.5
R <- R + t(R)
diag(R) <- 1

rep <- 50  # repeat times

##### LPM #####
library(LPM)

# function to generate data
generate_data <- function(M, K, beta0, alpha, R){
  
  Z <- matrix(rep(beta0, each = M), M, K) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- NULL
  
  for (k in 1:K){
    Pvalue_tmp <- runif(M)
    Pvalue_tmp[indexeta[, k]] <- rbeta(sum(indexeta[, k]), alpha[k], 1)
    
    Pvalue <- c(Pvalue, list(data.frame(SNP = seq(1, M), p = Pvalue_tmp)))
    
  }
  
  names(Pvalue) <- paste("P", seq(1, K), sep = "")
  
  return( list(Pvalue = Pvalue, eta = eta))
}

est_bLPM <- NULL
est_LPM <- NULL

for (i in 1:rep){
  data <- generate_data(M, K, beta0, alpha, R)
  Pvalue <- data$Pvalue

  fit <- bLPM(Pvalue, coreNum = 10)
  est <- c(est, list(fit))

  fitLPM <- LPM(fit)
  est_LPM <- c(est_LPM, list(fitLPM))
}

est_rho_LPM <- array(0, c(K, K, rep))
test_rho_LPM <- array(0, c(K, K, rep))
for(i in 1:rep){
  est_rho_LPM[, , i] <- est_LPM[[i]]$R
  rho_pvalue <- test_rho(est_blPM[[i]])
  test_rho_LPM[, , i] <- (rho_pvalue < 0.05/((K-1)*K/2))
}

# results to get Supplementary Figure S2a
test_rho_LPM <- apply(test_rho_LPM, c(1, 2), mean)
colnames(test_rho_LPM) <- paste("P", 1:8, sep = "")
rownames(test_rho_LPM) <- colnames(test_rho_LPM)


##### GPA #####
library(GPA)

# function to generate data
generate_data_GPA <- function(M, K, beta0, alpha, R){
  
  Z <- matrix(rep(beta0, each = M), M, K) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- matrix(runif(M*K), M, K)
  
  for (k in 1:K){
    Pvalue[indexeta[, k], k] <- rbeta(sum(indexeta[, k]), alpha[k], 1)
  }
  
  return( list(Pvalue = Pvalue, eta = eta))
}

test_rho_GPA <- array(0, c(K, K, rep))
for (k in 1:rep){
  data <- generate_data_GPA(M, K, beta0, alpha, R)
  Pvalue <- data$Pvalue

  for (i in 1:(K-1)){
    for (j in (i+1):K){
      fit <- GPA(Pvalue[, c(i, j)])
      fit.H0 <- GPA(Pvalue[, c(i, j)], pleiotropyH0 = TRUE)
      test <- pTest(fit, fit.H0)
      test_rho_GPA[i, j, k] <- (test$pvalue < 0.05/28)
    }
  }
}

# results to get Supplementary Figure S2b
test_rho_GPA <- apply(test_rho_GPA, c(1, 2), mean)
test_rho_GPA <- test_rho_GPA + t(test_rho_GPA)
diag(test_rho_GPA) <- 1
colnames(test_rho_GPA) <- paste("P", 1:8, sep = "")
rownames(test_rho_GPA) <- colnames(test_rho_GPA)

##### GGPA #####
library(GGPA)

# function to generate data
generate_data_GGPA <- function(M, K, beta0, alpha, R){
  
  Z <- matrix(rep(beta0, each = M), M, K) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- matrix(runif(M*K), M, K)
  
  for (k in 1:K){
    Pvalue[indexeta[, k], k] <- rbeta(sum(indexeta[, k]), alpha[k], 1)
  }
  
  return( list(Pvalue = Pvalue, eta = eta))
}

est_GGPA <- NULL

for (i in 1:rep){
  data <- generate_data_GGPA(M, K, beta0, alpha, R)

  fit_GGPA <- GGPA(data$Pvalue)

  est_GGPA <- c(est_GGPA, list(fit_GGPA))
}

# result to get Supplementary Figure S2c
plot(est_GGPA[[1]]) 
plot(est_GGPA[[2]])
