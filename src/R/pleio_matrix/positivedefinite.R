locisigma <- 0.1
covmean <- -0.5
pleio_covs <- rnorm(28, covmean, abs(covmean/20))

library(matrixcalc)

sigma <- matrix(c(locisigma, pleio_covs[0+1], pleio_covs[1+1], pleio_covs[2+1], pleio_covs[3+1], pleio_covs[4+1], pleio_covs[5+1], pleio_covs[6+1],
         pleio_covs[0+1], locisigma, pleio_covs[7+1], pleio_covs[8+1], pleio_covs[9+1], pleio_covs[10+1], pleio_covs[11+1], pleio_covs[12+1],
         pleio_covs[1+1], pleio_covs[7+1], locisigma, pleio_covs[13+1], pleio_covs[14+1], pleio_covs[15+1], pleio_covs[16+1], pleio_covs[17+1],
         pleio_covs[2+1], pleio_covs[8+1], pleio_covs[13+1], locisigma, pleio_covs[18+1], pleio_covs[19+1], pleio_covs[20+1], pleio_covs[21+1],
         pleio_covs[3+1], pleio_covs[9+1], pleio_covs[14+1], pleio_covs[18+1], locisigma, pleio_covs[22+1], pleio_covs[23+1], pleio_covs[24+1],
         pleio_covs[4+1], pleio_covs[10+1], pleio_covs[15+1], pleio_covs[19+1], pleio_covs[22+1], locisigma, pleio_covs[25+1], pleio_covs[26+1],
         pleio_covs[5+1], pleio_covs[11+1], pleio_covs[16+1], pleio_covs[20+1], pleio_covs[23+1], pleio_covs[25+1], locisigma, pleio_covs[27+1],
         pleio_covs[6+1], pleio_covs[12+1], pleio_covs[17+1], pleio_covs[21+1], pleio_covs[24+1], pleio_covs[26+1], pleio_covs[27+1], locisigma), nrow=8)

library(Matrix)

sigma_posdef <- nearPD(sigma, keepDiag = T)$mat
print(sigma_posdef)
write.matrix(sigma_posdef)

setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R")

{
  library(Matrix)
  s <- as.vector(t(read.csv("sig.csv", header = F)))
  sigma <- matrix(c(s), nrow=8)
  sigma_PD <- nearPD(sigma, keepDiag = T)$mat
  write.matrix(sigma_PD, file = "sig_PD.txt")
}

{
  library(Matrix)
  s <- as.vector(t(read.csv("sig.csv", header = F)))
  sigma <- matrix(c(s), nrow=8)
  sigma_PD <- nearPD(sigma, keepDiag = T)$mat
  write.matrix(sigma_PD, file = "PD_dx28d4.txt")
}
