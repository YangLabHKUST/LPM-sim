##### Performance for the hypothesis testing of annotation enrichment #####
# Vary alpha (0.2, 0.4, 0.6) and beta (−0.4, −0.3, −0.2, −0.1, 0.1, 0.2, 0.3, 0.4) to get Supplementary Figure S28
# Vary A.perc (0.001, 0.005, 0.01, 0.05, 0.1, 0.2) and beta (-0.3, -0.2, -0.1, 0.1, 0.2, 0.3) to get Supplementary Figures S29 and S30

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
beta    <- rep(0, K)   # coefficients of annotations
A.perc <- 0.2                        # the proportion the entries in X is 1
A         <- rep(0, M*D)             # the design matrix of annotation
indexA    <- sample(M*D, M*D*A.perc)
A[indexA] <- 1
A         <- matrix(A, M, D)

alpha <- 0.2   # parameter in the Beta distribution
rho <- 0       # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K) # correlation matrix for the traits

rep <- 500  # repeat times

pvalue_beta <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  fit <- bLPM(Pvalue, X = X)
  LPMfit <- LPM(fit)
  pvalue_beta[i] <- test_beta(Pvalue, X, 1, LPMfit)$p_value[2]
}

result <- sum(pvalue_beta < 0.05)/rep
