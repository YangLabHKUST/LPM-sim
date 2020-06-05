##### Parameter estimation #####
# Supplementary Figures S31-S49

library(MASS)
library(pbivnorm)
library(mvtnorm)

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

est_alpha <- NULL
est_beta <- NULL
est_R <- NULL

est_bLPM_alpha <- NULL
est_bLPM_beta <- NULL
est_bLPM_rho <- NULL

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A

  fit <- bLPM(Pvalue, X = X, coreNum = 10)
  fitLPM <- LPM(fit)
  
  est_alpha <- c(est_alpha, list(fitLPM$alpha))
  est_beta <- c(est_beta, list(fitLPM$beta))
  est_R <- c(est_R, list(fitLPM$R))
  
  est_bLPM_alpha <- c(est_bLPM_alpha, list(fit$alpha))
  est_bLPM_beta <- c(est_bLPM_beta, list(fit$beta))
  est_bLPM_rho <- c(est_bLPM_rho, list(fit$rho))
}

