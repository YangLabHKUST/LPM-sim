##### Performance the identification of risk SNPs for one trait (eight traits) #####
# Figure 2 in the main text and Supplementary Figure S18 and S19
# Change the target trait to get Supplementary Figures S8-S14

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

FDR1 <- matrix(0, rep, 9)
AUC1 <- matrix(0, rep, 9)
FDR14 <- matrix(0, rep, 7)
AUC14 <- matrix(0, rep, 7)
FDR124 <- numeric(rep)
AUC124 <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A

  fit <- bLPM(Pvalue, X = X, coreNum = 10)
  fitLPM <- LPM(fit)

  post <- post(Pvalue[1], X, 1, fitLPM)
  assoc1 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1[i, 1] <- comp_FDR(data$eta[, 1], assoc1$eta)
  AUC1[i, 1] <- comp_AUC(data$eta[, 1], post$posterior)
  
  for (k in 2:8){
    post <- post2(Pvalue[c(1, i)], X, c(1, i), fitLPM)
    assoc2 <- assoc(post, FDRset = 0.1, fdrControl = "global")
    FDR1[i, k] <- comp_FDR(data$eta[, 1], assoc2$eta.marginal1)
    AUC1[i, k] <- comp_AUC(data$eta[, 1], post$post.marginal1)
  }
  
  post <- post3(Pvalue[1:3], X, c(1, 2, 3), fitLPM)
  assoc3 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1[i, 9] <- comp_FDR(data$eta[, 1], assoc3$eta.marginal1)
  AUC1[i, 9] <- comp_AUC(data$eta[, 1], post$post.marginal1)
  
  post <- post(Pvalue[c(1, 4)], X, c(1, 4), fitLPM)
  assoc2 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR14[i, 1] <- comp_FDR(((data$eta[, 1] + data$eta[, 4]) == 2), assoc2$eta.joint)
  AUC14[i, 1] <- comp_AUC(((data$eta[, 1] + data$eta[, 4]) == 2), post$post.joint)
  
  post <- post(Pvalue[c(1, 2, 4)], X, c(1, 2, 4), fitLPM)
  assoc3 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR14[i, 2] <- comp_FDR(((data$eta[, 1] + data$eta[, 4]) == 2), assoc3$eta.marginal13)
  AUC14[i, 2] <- comp_AUC(((data$eta[, 1] + data$eta[, 4]) == 2), post$post.marginal13)
  FDR124[i] <- comp_FDR(((data$eta[, 1] + data$eta[, 2] + data$eta[, 4]) == 3), assoc3$eta.joint)
  AUC124[i] <- comp_AUC(((data$eta[, 1] + data$eta[, 2] + data$eta[, 4]) == 3), post$post.joint)
  
  id <- 3
  for (k in c(3, 5:8)){
    post <- post(Pvalue[c(1, 4, k)], X, c(1, 4, k), fitLPM)
    assoc3 <- assoc(post, FDRset = 0.1, fdrControl = "global")
    FDR14[i, id] <- comp_FDR(((data$eta[, 1] + data$eta[, 4]) == 2), assoc3$eta.marginal12)
    AUC14[i, id] <- comp_AUC(((data$eta[, 1] + data$eta[, 4]) == 2), post$post.marginal12)
    id <- id + 1
  }
  
}

##### GPA #####
library(GPA)

# function to generate data
generate_data_GPA <- function(M, K, D, A, beta, alpha, R){
  
  Z <- cbind(rep(1, M), A) %*% t(beta) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- matrix(0, M, K)
  
  for (k in 1:K){
    Pvalue[, k] <- runif(M)
    Pvalue[indexeta[, k], k] <- rbeta(sum(indexeta[, k]), alpha[k], 1)
  }
  
  return( list(Pvalue = Pvalue, A = A, beta = beta, eta = eta))
}

FDR1_GPA <- matrix(0, rep, 9)
AUC1_GPA <- matrix(0, rep, 9)
FDR14_GPA <- matrix(0, rep, 7)
AUC14_GPA <- matrix(0, rep, 7)
FDR124_GPA <- numeric(rep)
AUC124_GPA <- numeric(rep)

for (i in 1:rep){
  data <- generate_data_GPA(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A

  fit <- GPA(Pvalue[, 1], X)
  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global")
  FDR1_GPA[i, 1] <- comp_FDR(data$eta[, 1], assoc1)
  AUC1_GPA[i, 1] <- comp_AUC(data$eta[, 1], fdr(fit))
  
  for (k in 2:8){
    fit <- GPA(Pvalue[, c(1, k)], X)
    assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "1*")
    FDR1_GPA[i, k] <- comp_FDR(data$eta[, 1], assoc1)
    AUC1_GPA[i, k] <- comp_AUC(data$eta[, 1], fdr(fit, pattern = "1*"))
  }
  
  fit <- GPA(Pvalue[, c(1, 2, 3)], X)
  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "1**")
  FDR1_GPA[i, 9] <- comp_FDR(data$eta[, 1], assoc1)
  AUC1_GPA[i, 9] <- comp_AUC(data$eta[, 1], fdr(fit, pattern = "1**"))
  
  fit <- GPA(Pvalue[, c(1, 4)], X)
  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "11")
  FDR14_GPA[i, 1] <- comp_FDR(((data$eta[, 1] + data$eta[, 4]) == 2), assoc1)
  AUC14_GPA[i, 1] <- comp_AUC(((data$eta[, 1] + data$eta[, 4]) == 2), fdr(fit, pattern = "11"))
  
  id <- 1
  for (k in c(2, 3, 5:8)){
    id <- id + 1
    fit <- GPA(Pvalue[, c(1, 4, k)], X)
    assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "11*")
    FDR14_GPA[i, id] <- comp_FDR(((data$eta[, 1] + data$eta[, 4]) == 2), assoc1)
    AUC14_GPA[i, id] <- comp_AUC(((data$eta[, 1] + data$eta[, 4]) == 2), fdr(fit, pattern = "11*"))
  }
  
  fit <- GPA(Pvalue[, c(1, 2, 4)], X)
  assoc1 <- assoc(fit, FDR = 0.1, fdrControl = "global", pattern = "111")
  FDR124_GPA[i] <- comp_FDR(((data$eta[, 1] + data$eta[, 2] + data$eta[, 4]) == 3), assoc1)
  AUC124_GPA[i] <- comp_AUC(((data$eta[, 1] + data$eta[, 2] + data$eta[, 4]) == 3), fdr(fit, pattern = "111"))
  
}


