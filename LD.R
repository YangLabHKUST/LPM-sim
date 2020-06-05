##### Simulation study for evaluating the LD effects on LPM #####
# Supplementary Figure S49

####### prepare ##########
# download data from WTCCC https://www.wtccc.org.uk/info/access_to_data_samples.html

# quality control
system("plink-1.07-mac-intel/plink --bfile ctrl --geno 0.01 --mind 0.05 --hwe 0.001 --maf 0.05 --make-bed --out ctrlqc")

# use only SNPs in CHR 1
system("plink-1.07-mac-intel/plink --bfile ctrlqc --chr 1 --make-bed --out ctrlqc_chr1")

# keep the SNPs with BP
SNP <- read.table(file = "ctrlqc_chr1.bim")
ref <- read.table(file = "wtccc.variant_function")
ref_chr1 <- ref[which(ref$V3 == "chr1"), ]
save(ref_chr1, file = "ref_chr1.RData")
keep.SNP <- as.vector(SNP$V2[which(SNP$V2 %in% ref_chr1$V8)])
write.table(keep.SNP, file = "keep.SNP", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("plink-1.07-mac-intel/plink --bfile ctrlqc_chr1 --extract keep.SNP --make-bed --out ctrl_qc_chr1")

# keep overlaped SNPs with Annotation
SNP <- read.table(file = "ctrl_qc_chr1.bim")
load("region9tissue127.RData") # Annotation data

SNP_share <- as.character(SNP$V2[SNP$V2 %in% region9tissue127$SNP])
rownames(region9tissue127) <- region9tissue127$SNP
Annot <- region9tissue127[SNP_share, ]
save(Annot, file = "Annot.RData")

write.table(SNP_share, file = "SNP_share", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("plink-1.07-mac-intel/plink --bfile ./used/ctrl_qc_chr1 --extract SNP_share --make-bed --out ctrl_qc_chr1")

# seperate individuals
ind <- read.table(file = "ctrl_qc_chr1.fam")
set.seed(1)
index <- sample(nrow(ind), nrow(ind)/2)
ind <- ind[index, 1:2]
write.table(ind, file = "ind_filter.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("plink-1.07-mac-intel/plink --bfile ctrl_qc_chr1 --keep ind_filter.txt --make-bed --out ctrl_qc_chr1_1")
system("plink-1.07-mac-intel/plink --bfile ctrl_qc_chr1 --remove ind_filter.txt --make-bed --out ctrl_qc_chr1_2")

# set causal SNPs
library(MASS)
load("Annot.RData")
Annot <- Annot[, c(1:10)]
D <- ncol(Annot) - 1
M <- nrow(Annot)
r <- 0.25                        # the relative signal strength between annotated part and un-annotated part
A <- as.matrix(Annot[, -1])

beta0 <- -3.54           # intercept of the probit model
set.seed(1)
beta    <- rnorm(D)  # coefficients of annotations
sigmae2 <- var(A %*% as.matrix(beta))/r
beta    <- beta/sqrt(diag(sigmae2))
beta    <- c(beta0, beta)

set.seed(1)
Z <- cbind(rep(1, M), A) %*% as.matrix(beta) + rnorm(M)
indexeta <- (Z > 0)

causal.SNP1 <- rownames(Annot)[which(indexeta[, 1] == 1)]
write.table(causal.SNP1, file = "causal1.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)

beta0 <- -3.501           # intercept of the probit model
set.seed(2)
beta    <- rnorm(D)  # coefficients of annotations
sigmae2 <- var(A %*% as.matrix(beta))/r
beta    <- beta/sqrt(diag(sigmae2))
beta    <- c(beta0, beta)
set.seed(2)
Z <- cbind(rep(1, M), A) %*% as.matrix(beta) + rnorm(M)
indexeta <- (Z > 0)

causal.SNP2.0 <- rownames(Annot)[which(indexeta[, 1] == 1)] # 100% shared risk SNPs
write.table(causal.SNP2.0, file = "causal2.0_new.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)

causal.SNP2.05 <- causal.SNP2.0 # 50% shared risk SNPs
causal.SNP2.05[sample(length(causal.SNP2.0), round(length(causal.SNP2.0)/2))] <- causal.SNP1[sample(length(causal.SNP2.0), round(length(causal.SNP2.0)/2))]
write.table(causal.SNP2.05, file = "causal2.05_new.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)
causal.SNP2.1 <- causal.SNP1 # 100% shared risk SNPs
write.table(causal.SNP2.1, file = "causal2.1.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)

# simulate GWAS data
set.seed(1)
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_1  --simu-qt  --simu-causal-loci causal1.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl1")
set.seed(1)
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_2  --simu-qt  --simu-causal-loci causal2.0.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl2.0")
set.seed(1)
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_2  --simu-qt  --simu-causal-loci causal2.05.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl2.05")
set.seed(1)
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_2  --simu-qt  --simu-causal-loci causal2.1.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl2.1")

# copy ctrl_qc_chr1_1 and rename it as tmp1
# copy ctrl_qc_chr1_2 and rename it as tmp2

library(LPM)
tmp1 <- read.table(file = "tmp1.fam")
tmp1_bim <- read.table("tmp1.bim")
sim1 <- read.table(file = "sim_ctrl1.phen")
SNP <- read.table(file = "tmp1.bim")
tmp2 <- read.table(file = "tmp2.fam")
tmp2_bim <- read.table("tmp2.bim")
sim2.0 <- read.table(file = "sim_ctrl2.0_new.phen")
sim2.05 <- read.table(file = "sim_ctrl2.05_new.phen")
sim2.1 <- read.table(file = "sim_ctrl2.1.phen")

load("Annot.RData")
A <- Annot[, 1:10]

M <- nrow(A)

est.joint.0 <- matrix(0, 50, M)
est.joint.05 <- matrix(0, 50, M)
est.joint.1 <- matrix(0, 50, M)

for (i in 1:50){
  # no pleiotropy
  tmp1$V6 <- sim1[, i+2]
  write.table(tmp1, file = "tmp1.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp1 --assoc  --out sim_pvalue1")
  pvalue1 <- read.table(file = "sim_pvalue1.qassoc", header = TRUE)
  tmp2$V6 <- sim2.0[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.0")
  pvalue2.0 <- read.table(file = "sim_pvalue2.0.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.0$SNP, p = pvalue2.0$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.joint.0[i, ] <- assoc.SNP$eta.marginal1
  # P2 50% shared risk SNPs
  tmp2$V6 <- sim2.05[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.05")
  pvalue2.05 <- read.table(file = "sim_pvalue2.05.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.05$SNP, p = pvalue2.05$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.joint.05[i, ] <- assoc.SNP$eta.marginal1
  # P2 100% shared risk SNPs
  tmp2$V6 <- sim2.1[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.1")
  pvalue2.1 <- read.table(file = "sim_pvalue2.1.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.1$SNP, p = pvalue2.1$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.joint.1[i, ] <- assoc.SNP$eta.marginal1
}

load("ref_chr1.RData")
SNP <- read.table(file = "tmp1.bim")

FDR.summary <- NULL

pleiotropy <- c("no pleiotropy", "50% shared risk SNPs", "100% shared risk SNPs")
plei <- c("0", "05", "1")

distance <- seq(1e5, 1e6, 1e5)

for(j in 1:length(plei)){
  causal.SNP <- read.table(file = paste("causal1.snplist", sep = ""))
  causal.bp <- ref_chr1$V4[which((ref_chr1$V8 %in% causal.SNP$V1) == 1)]
  
  est <- get(paste("est.joint.", plei[j], sep = ""))
  FDR <- matrix(0, 50, 10)
  
  for (k in 1:length(distance)){
    bp.range <- NULL
    for (l in 1:length(causal.bp)){
      bp.range <- c(bp.range, seq(causal.bp[l] - distance[k], causal.bp[l] + distance[k]))
    }
    bp.range <- unique(bp.range)
    SNP.range <- ref_chr1$V8[which((ref_chr1$V4 %in% bp.range) == 1)]
    SNP.range <- SNP.range[which((SNP.range %in% SNP$V2) == 1)]
    for (l in 1:50){
      TP <- sum(SNP$V2[which(est[l, ] == 1)] %in% SNP.range)
      
      if (sum(est[l, ]) == 0)
        FDR[l, k] <- 0
      else
        FDR[l, k] <- (sum(est[l, ]) - TP)/sum(est[l, ])
    }
  }
  
  FDR.summary$pleiotropy <- c(FDR.summary$pleiotropy, rep(pleiotropy[j], 50*10))
  FDR.summary$distance <- c(FDR.summary$distance, rep(distance, each = 50))
  FDR.summary$FDR <- c(FDR.summary$FDR, as.vector(FDR))
  
}


FDR.summary <- as.data.frame(FDR.summary)
FDR.summary$distance <- as.factor(FDR.summary$distance)
levels(FDR.summary$distance) <- c(paste(seq(100,1000,100), sep = ""))
FDR.summary$pleiotropy <- factor(FDR.summary$pleiotropy, levels = levels(FDR.summary$pleiotropy)[c(3, 2, 1)])
