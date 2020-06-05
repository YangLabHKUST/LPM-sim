##### Simulations if P-values are not from beta distribution #####
# Vary dist ("near_normal", "skew", "big_normal") and rho (0, 0.2, 0.4, 0.6) 
# to get Supplementary Figures S45 and S46

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

spiky <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 0.25)
  if(sum(r <= 0.2) != 0)
    z[which(r <= 0.2)] <- rnorm(sum(r <= 0.2), 0, 0.5)
  if(sum(r > 0.2 & r <= 0.4) != 0)
    z[which(r > 0.2 & r <= 0.4)] <- rnorm(sum(r > 0.2 & r <= 0.4), 0, 1)
  if(sum(r > 0.8) != 0)
    z[which(r > 0.8)] <- rnorm(sum(r > 0.8), 0, 2)
  
  return(z)
}

near_normal <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 1)
  if(sum(r <= 1/3) != 0)
    z[which(r <= 1/3)] <- rnorm(sum(r <= 1/3), 0, 2)
  
  return(z)
}

skew <- function(N){
  r <- runif(N)
  z <- rnorm(N, -2, 2)
  if(sum(r <= 0.25) != 0)
    z[which(r <= 0.25)] <- rnorm(sum(r <= 0.25), -1, 1.5)
  if(sum(r > 0.25 & r <= (0.25+1/3)) != 0)
    z[which(r > 0.25 & r <= (0.25+1/3))] <- rnorm(sum(r > 0.25 & r <= (0.25+1/3)), 0, 1)
  if(sum(r > 5/6) != 0)
    z[which(r > 5/6)] <- rnorm(sum(r > 5/6), 1, 1)
  
  return(z)
}

big_normal <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 4)
  
  return(z)
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

library(LPM)

# function to generate data
generate_data <- function(M, K, D, A, beta, R, dist){
  
  Z <- cbind(rep(1, M), A) %*% t(beta) + mvrnorm(M, rep(0, K), R)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, K)
  
  Pvalue <- NULL
  
  fz <- match.fun(dist)
  
  for (k in 1:K){
    z <- fz(sum(indexeta[, k]))
    Pvalue_tmp <- runif(M)
    Pvalue_tmp[indexeta[, k]] <- pnorm(abs(z), lower.tail = FALSE)*2
    
    Pvalue <- c(Pvalue, list(data.frame(SNP = seq(1, M), p = Pvalue_tmp)))
    
  }
  
  names(Pvalue) <- paste("P", seq(1, K), sep = "")
  
  A <- data.frame(SNP=seq(1,M), A)
  
  return( list(Pvalue = Pvalue, A = A, beta = beta, eta = eta))
}

# compute type I error
rho <- 0        # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K)
dist <- "spiky" # distribution
rep <- 500      # repeat times
pvalue_rho <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, R, dist)
  Pvalue <- data$Pvalue
  X      <- data$A

  fit <- bLPM(Pvalue, X = X)
  pvalue_rho[i] <- test_rho(fit)
}

TypeIerror <- sum(pvalue_rho < 0.05)/rep

# compute FDR
rho <- 0        # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K)
dist <- "spiky" # distribution
rep <- 50       # repeat times

FDR1.sep <- numeric(rep)
FDR2.sep <- numeric(rep)
FDR1.joint <- numeric(rep)
FDR2.joint <- numeric(rep)
FDR12 <- numeric(rep)

for (i in 1:rep){
  data <- generate_data(M, K, D, A, beta, R, dist)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  fit <- bLPM(Pvalue, X = X)
  fitLPM <- getLPMest(fit)
  
  post <- post(Pvalue[1], X, 1, fitLPM)
  assoc1 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1.sep[i] <- comp_FDR(data$eta[, 1], assoc1$eta)
  
  post <- post(Pvalue[2], X, 2, fitLPM)
  assoc1 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR2.sep[i] <- comp_FDR(data$eta[, 2], assoc1$eta)
  
  post <- post(Pvalue[c(1, 2)], X, c(1, 2), fitLPM)
  assoc2 <- assoc(post, FDRset = 0.1, fdrControl = "global")
  FDR1.joint[i] <- comp_FDR(data$eta[, 1], assoc2$eta.marginal1)
  FDR2.joint[i] <- comp_FDR(data$eta[, 2], assoc2$eta.marginal2)
  FDR12[i] <- comp_FDR(((data$eta[, 1] + data$eta[, 2]) == 2), assoc2$eta.joint)
}
