##### Performance in characterizing the correlation among the traits (two traits) #####
# Vary alpha (0.2, 0.4, 0.6), r (0.25, 1, 4) and rho (0, 0.05, 0.1, 0.15, 0.2, 0.25) 
# to get Supplementary Figure S7

library(MASS)
library(LPM)
library(pbivnorm)
library(mvtnorm)

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

alpha <- 0.2   # parameter in the Beta distribution
rho <- 0       # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K) # correlation matrix for the traits

rep <- 500  # repeat times

pvalue_rho <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  fit <- bLPM(Pvalue, X = X)
  pvalue_rho[i] <- test_rho(fit)
}

result <- sum(pvalue_rho < 0.05)/rep
