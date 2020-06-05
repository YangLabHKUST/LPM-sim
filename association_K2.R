##### Performance the identification of risk SNPs for one trait (two traits) #####
# Vary alpha2 (0.2, 0.4, 0.6) and rho (0, 0.2, 0.4, 0.6) to get Supplementary Figure S15

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

K <- 2                 # No. of traits
M <- 100000            # No. of SNPs
D <- 5                 # No. of annotations
beta0 <- -1            # intercept of the probit model
beta0 <- rep(beta0, K)
set.seed(1)
beta    <- matrix(rnorm(K*D), K, D)  # coefficients of annotations
A.perc <- 0.2                        # the proportion the entries in X is 1
A         <- rep(0, M*D)             # the design matrix of annotation
indexA    <- sample(M*D, M*D*A.perc)
A[indexA] <- 1
A         <- matrix(A, M, D)
r <- 1                               # the relative signal strengh between annotated part and un-annotated part
sigmae2 <- var(A %*% t(beta))/r
beta    <- beta/sqrt(diag(sigmae2))
beta    <- cbind(as.matrix(beta0), beta)

alpha1 <- 0.2  # parameter in the Beta distribution for the first trait
alpha2 <- 0.2  # parameter in the Beta distribution for the second trait
rho <- 0       # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K) # correlation matrix for the traits

rep <- 50  # repeat times

##### LPM #####
library(LPM)

# function to generate data
generate_data <- function(M, K, D, A, beta, alpha, R){
  
  Z <- cbind(rep(1, M), A) %*% t(beta) + mvrnorm(M, rep(0, K), R)

  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)

  Pvalue <- NULL

  for (k in 1:K){
    Pvalue_tmp <- runif(M)
    Pvalue_tmp[indexeta[, k]] <- rbeta(sum(indexeta[, k]), alpha[k], 1)

    Pvalue <- c(Pvalue, list(data.frame(SNP = seq(1, M), p = Pvalue_tmp)))

  }

  names(Pvalue) <- paste("P", seq(1, K), sep = "")

  A <- data.frame(SNP=seq(1,M), A)

  return( list(Pvalue = Pvalue, A = A, beta = beta, eta = eta))
}

FDR <- numeric(rep)
AUC <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A

  fit <- bLPM(Pvalue, X = X)
  fitLPM <- LPM(fit)

  post <- post(Pvalue, X, 1:2, fitLPM)
  assoc <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR[i] <- comp_FDR(data$eta[, 1], assoc$eta.marginal1)
  AUC[i] <- comp_AUC(data$eta[, 1], post$post.marginal1)
}
