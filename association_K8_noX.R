##### Performance the identification of risk SNPs for one trait when annotations have no role (eight traits) #####
# Change the target trait to get Figures S20-S27 in Supplementary Document

library(MASS)
library(pbivnorm)
library(mvtnorm)
library(pROC)

# function to compute FDR
comp_FDR <- function(true, est){
  
  t <- table(true, est)
  if (sum(est)==0){
    FDR.fit <- 0
  }
  else if (sum(est)==length(est)){
    FDR.fit <- t[1]/(t[1]+t[2])
  }
  else{
    FDR.fit <- t[1,2]/(t[1,2]+t[2,2])
  }
  
  return(FDR.fit)
}
# function to compute AUC
comp_AUC <- function(true, post){
  fdr <- 1 - post
  AUC <- as.numeric(roc(true, fdr)$auc)
  return(AUC)
}

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

FDR1 <- matrix(0, rep, 9)
AUC1 <- matrix(0, rep, 9)

for (i in 1:rep){
  data <- generate_data(M, K, beta0, alpha, R)
  Pvalue <- data$Pvalue

  fit <- bLPM(Pvalue, X = NULL, coreNum = 10)
  fitLPM <- LPM(fit)

  post <- post(Pvalue[1], X = NULL, 1, fitLPM)
  assoc1 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1[i, 1] <- comp_FDR(data$eta[, 1], assoc1$eta)
  AUC1[i, 1] <- comp_AUC(data$eta[, 1], post$posterior)
  
  for (k in 2:8){
    post <- post2(Pvalue[c(1, k)], X = NULL, c(1, k), fitLPM)
    assoc2 <- assoc(post, FDRset = 0.1, fdrControl = "global")
    FDR1[i, k] <- comp_FDR(data$eta[, 1], assoc2$eta.marginal1)
    AUC1[i, k] <- comp_AUC(data$eta[, 1], post$post.marginal1)
  }
  
  post <- post3(Pvalue[1:3], X = NULL, c(1, 2, 3), fitLPM)
  assoc3 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1[i, 9] <- comp_FDR(data$eta[, 1], assoc3$eta.marginal1)
  AUC1[i, 9] <- comp_AUC(data$eta[, 1], post$post.marginal1)
}

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

FDR1_GPA <- matrix(0, rep, 9)
AUC1_GPA <- matrix(0, rep, 9)

for (i in 1:rep){
  data <- generate_data_GPA(M, K, beta0, alpha, R)
  Pvalue <- data$Pvalue

  fit <- GPA(Pvalue[, 1])
  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global")
  FDR1_GPA[i, 1] <- comp_FDR(data$eta[, 1], assoc1)
  AUC1_GPA[i, 1] <- comp_AUC(data$eta[, 1], fdr(fit))
  
  for (k in 2:8){
    fit <- GPA(Pvalue[, c(1, k)])
    assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "1*")
    FDR1_GPA[i, k] <- comp_FDR(data$eta[, 1], assoc1)
    AUC1_GPA[i, k] <- comp_AUC(data$eta[, 1], fdr(fit, pattern = "1*"))
  }
  
  fit <- GPA(Pvalue[, c(1, 2, 3)])
  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "1**")
  FDR1_GPA[i, 9] <- comp_FDR(data$eta[, 1], assoc1)
  AUC1_GPA[i, 9] <- comp_AUC(data$eta[, 1], fdr(fit, pattern = "1**"))
  
}


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

FDR_GGPA <- matrix(0, rep, 8)
AUC_GGPA <- matrix(0, rep, 8)

for (i in 1:rep){
  data <- generate_data_GGPA(M, K, beta0, alpha, R)

  fit_GGPA <- GGPA(data$Pvalue)

  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global")
  for (k in 1:8){
    FDR[i, k] <- comp_FDR(data$eta[, k], assoc1[, k])
    power[i, k] <- comp_power(data$eta[, k], assoc1[, k])
  }
}

