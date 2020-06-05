##### Performance when the proportion of risk SNPs is extremely small #####
# Vary alpha (0.2, 0.4, 0.6), pi1 (0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2) and rho (0, 0.2, 0.4, 0.6) 
# to get Supplementary Figures S47 and S48

library(MASS)
library(LPM)
library(pbivnorm)
library(mvtnorm)

# function to generate data
generate_data <- function(M, K, beta0, alpha, R){
  
  Z <- matrix(rep(beta0, each = M), M, K) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- NULL
  
  for (k in 1:K){
    Pvalue_tmp <- runif(M)
    Pvalue_tmp[indexeta[, k]] <- rbeta(sum(indexeta[, k]), alpha, 1)
    
    Pvalue <- c(Pvalue, list(data.frame(SNP = seq(1, M), p = Pvalue_tmp)))
    
  }
  
  names(Pvalue) <- paste("P", seq(1, K), sep = "")
  
  return(list(Pvalue = Pvalue, eta = eta))
}

K <- 2                   # No. of traits
M <- 100000              # No. of SNPs
pi1 <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2) # proportion of risk SNPs
beta0 <- -qnorm(1 - pi1) # intercept of the probit model

alpha <- c(0.2, 0.4, 0.6)   # parameter in the Beta distribution
rho <- c(0, 0.2, 0.4, 0.6)  # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K) # correlation matrix for the traits

rep <- 50  # repeat times

est_pi1 <- numeric(rep)
est_rho <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, beta0, rep(alpha, K), R)
  Pvalue <- data$Pvalue
  
  fit <- bLPM(Pvalue, X = NULL)
  est_pi1[i] <- pnorm(fit$beta0[1])
  est_rho[i] <- fit$rho
}

