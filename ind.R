##### Simulations when P-values are obtained from individual-level data #####
# Vary N2 (3000, 5000, 10000), h2 (0.3, 0.5, 0.8) and rho (0, 0.2, 0.4, 0.6, 0.8) 
# to get Supplementary Figures S42-S44

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

K <- 2                # No. of traits
M <- 20000            # No. of SNPs
N1 <- 10000           # No. of individuals for the first GWAS
N2 <- 3000            # No. of individuals for the first GWAS
N <- c(N1, N2)   
D <- 5                # No. of annotations
beta0 <- -3           # intercept of the probit model
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
f <- runif(M, 0.05, 0.5)             # minor allele frequency

library(LPM)

# function to generate data
generate_data_ind <- function(M, K, D, N, A, beta, h2, R, f){
  
  Z <- cbind(rep(1, M), A) %*% t(beta) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- NULL
  
  for (k in 1:K){
    # genotype data
    X_tmp <- matrix(rnorm(N[k]*M), N[k], M)
    X <- matrix(1, N[k], M)
    X[t(t(X_tmp) - quantile(X_tmp, f^2)) < 0] <- 2
    X[t(t(X_tmp) - quantile(X_tmp, 1-(1-f)^2)) > 0] <- 0
    
    # effect size
    beta_SNP <-  numeric(M)
    beta_SNP[indexeta[, k]] <- rnorm(sum(indexeta[, k]), 0, 1)
    
    # environment effect
    e <- rnorm(N, 0, sqrt((1/h2-1)*var(X%*%beta_SNP)))
    
    # phenotype
    y <- X%*%beta_SNP + e
    
    # p-value
    Pvalue_tmp <- apply(X, 2, function(X.col) summary(lm(y ~ X.col))$coefficients[2,4])
    
    Pvalue <- c(Pvalue, list(data.frame(SNP = seq(1, M), p = Pvalue_tmp)))
    
  }
  
  names(Pvalue) <- paste("P", seq(1, K), sep = "")
  
  A <- data.frame(SNP=seq(1,M), A)
  
  return( list(Pvalue = Pvalue, A = A, beta = beta, eta = eta))
}

# compute type I error
rho <- 0    # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K)
h2 <- 0.3   # heritability
rep <- 500  # repeat times

pvalue_rho <- numeric(rep)

for (i in 1:rep){
  data <- generate_data_ind(M, K, D, N, A, beta, h2, R, f)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  fit <- bLPM(Pvalue, X = X)
  pvalue_rho[i] <- test_rho(fit)[1, 2]
}

TypeIerror <- sum(pvalue_rho < 0.05)/rep

# estimate alpha
rho <- 0  # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K)
h2 <- 0.3 # heritability
rep <- 50 # repeat times

est_alpha <- numeric(rep)

for (i in 1:rep){
  data <- generate_data_ind(M, K, D, N, A, beta, h2, R, f)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  fit <- bLPM(Pvalue, X = X)
  est_alpha[i] <- fit$alpha[2, 1]
}

# compute FDR
rho <- 0  # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K)
h2 <- 0.3 # heritability
rep <- 50 # repeat times

FDR1 <- numeric(rep)
FDR2 <- numeric(rep)
FDR12 <- numeric(rep)

for (i in 1:rep){
  data <- generate_data_ind(M, K, D, N, A, beta, h2, R, f)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  fit <- bLPM(Pvalue, X = X)
  fitLPM <- LPM(fit)
  
  post <- post(Pvalue[c(1, 2)], X, c(1, 2), fitLPM)
  assoc2 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1[i] <- comp_FDR(data$eta[, 1], assoc2$eta.marginal1)
  FDR2[i] <- comp_FDR(data$eta[, 2], assoc2$eta.marginal2)
}
