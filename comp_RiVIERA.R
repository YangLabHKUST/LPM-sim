##### Comparison with RiVIERA #####
# Supplementary Figures S16 and S17

library(MASS)
library(pbivnorm)
library(mvtnorm)

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

alpha <- c(0.2, 0.2)  # parameter in the Beta distribution
rho <- 0.6            # correlation between the two traits
R <- matrix(c(1, rho, rho, 1), K, K) # correlation matrix for the traits

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

eta_all <- NULL
post_all <- NULL
ppa_all <- NULL
pvalue_all <- NULL

for (i in 1:rep){
  
  data <- generate_data(M, K, D, A, beta, alpha, R)
  Pvalue <- data$Pvalue
  X      <- data$A
  
  eta_all <- c(eta_all, list(data$eta))

  # LPM
  fit <- bLPM(Pvalue, X = X)
  fitLPM <- LPM(fit)
  post <- post(Pvalue[c(1, 2)], X, c(1, 2), fitLPM)
  post_all <- c(post_all, list(post))

  # RiVIERA
  gwasPval <- cbind(Pvalue[[1]]$p, Pvalue[[2]]$p)
  annotMat <- as.matrix(X[, -1])
  locus_cnt <- rep(1000, M/1000)
  locus_index <- make_locus_index(locus_cnt = locus_cnt)
  
  nsteps <- 100
  max_epoch <- 1e3
  step <- 1e-3
  burnFrac <- 0.2
  set.seed(100)
  ensemble_fit <- rivieraBeta(gwasPval=gwasPval,
                              annMat=annotMat,
                              positiveAnnotConstraint=TRUE,
                              locus_index=locus_index,
                              max_epoch=max_epoch,
                              step=step, nsteps=nsteps,
                              printfreq=10, verbose=FALSE)
  
  burnin <- round(burnFrac * slot(ensemble_fit, "fit_info")$ensembleSize)
  
  riviera_ppa <- finemap(ensemble=ensemble_fit,
                         gwasPval=gwasPval,
                         annMat = annotMat,
                         locus_index=matrix(c(1, 0, 99999), 1, 3),
                         burnin=burnin)
  ppa_all <- c(ppa_all, list(riviera_ppa))

  # pvalue
  pvalue_all <- c(pvalue_all, list(Pvalue))
}    

top <- c(10000, 15000, 18000, 20000)
top_LPM <- matrix(0, rep, length(top))
for(i in 1:rep){
  for(j in 1:length(top)){
    post_tmp <- post_all[[i]]$post.marginal1
    est_tmp <- (rank(post_tmp) > (length(post_tmp)-top[j]))
    eta_tmp <- eta_all[[i]][, 1]
    top_LPM[i, j] <- sum(((est_tmp + eta_tmp)==2))/sum(eta_tmp)
  }
}

top_riviera <- matrix(0, rep, length(top))
for(i in 1:rep){
  for(j in 1:length(top)){
    post_tmp <- ppa_all[[i]][, 1]
    est_tmp <- (rank(post_tmp) > (length(post_tmp)-top[j]))
    eta_tmp <- eta_all[[i]][, 1]
    top_riviera[i, j] <- sum(((est_tmp + eta_tmp)==2))/sum(eta_tmp)
  }
}

top_pvalue <- matrix(0, rep, length(top))
for(i in 1:rep){
  for(j in 1:length(top)){
    pvalue_tmp <- pvalue_all[[i]]$P1$p
    est_tmp <- (rank(-pvalue_tmp) > (length(pvalue_tmp)-top[j]))
    eta_tmp <- eta_p_all[[i]][, 1]
    top_pvalue[i, j] <- sum(((est_tmp + eta_tmp)==2))/sum(eta_tmp)
  }
}

eta1 <- matrix(0, rep, M)
comp_roc_LPM <- matrix(0, rep, M)
comp_roc_riviera <- matrix(0, rep, M)
for (i in 1:rep){
  eta1[i, ] <- eta_all[[i]][, 1]
  comp_roc_LPM[i, ] <- post_all[[i]]$post.marginal1
  comp_roc_riviera[i, ] <- ppa_all[[i]][, 1]
}

library(ROCR)
x <- prediction(as.list(as.data.frame(t(rbind(comp_roc_LPM, comp_roc_riviera)))), as.list(as.data.frame(t(rbind(eta1, eta1)))))
ROC <- performance(x, "tpr", "fpr")
