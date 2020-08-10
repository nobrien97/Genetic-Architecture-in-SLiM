# Pilot data analysis: Take G matrices, eigenanalysis, tensors? 

setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs")

# Basic standard error function for use with data manipulation

std.error <- function(x){
  n <- length(x)
  sd <- sd(x)
  sd/sqrt(n)
}




# Null data set

d_null <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/pilot_runs/Output/out_8T_null_means.csv", header = F)

# Name columns

names(d_null)[1:7]<- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu", "nloci")


# pleiocov terms                
names(d_null)[8:35] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))


names(d_null)[36] <- "Ne"

names(d_null)[37:44] <- paste0("mean", 0:7)

names(d_null)[45:52] <- paste0("var", 0:7)

names(d_null)[53:80] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null)[81:108] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null)[109:118] <- paste0("H_chr", 0:9)

# Organise data for G matrix construction

library(dplyr)

d_sbstfinal_null <- as_tibble(d_null[d_null$gen == 150000,])

d_null_mat <- select(d_sbstfinal_null, gen, seed, modelindex, paste0("var", 0:7), c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7)))

d_null_mat <- d_null_mat %>%
  group_by(gen, modelindex) %>%
  summarise_all(list(groupmean = mean, se = std.error))


# Data coercion to matrix function: expects a single row of 39 variables, gen, seed, modelindex, 8 variances, 28 covariances

dat_to_mat <- function(dat) {
  dat <- as.vector(t(dat))
  vars <- dat[4:11]
  covs <- dat[12:39]
  M <- matrix(data = c(
    vars[1], covs[1:7],
    covs[1], vars[2], covs[8:13],
    covs[c(2, 8)], vars[3], covs[14:18],
    covs[c(3, 9, 14)], vars[4], covs[19:22],
    covs[c(4, 10, 15, 19)], vars[5], covs[23:25],
    covs[c(5, 11, 16, 20, 23)], vars[6], covs[26:27],
    covs[c(6, 12, 17, 21, 24, 26)], vars[7], covs[28],
    covs[c(7, 13, 18, 22, 25, 27, 28)], vars[8]), nrow = 8, byrow = T
    )
  M
}

# Split data frame up into single lines, do dat_to_mat on each, store as an array

mat_gen <- function(dat) {
  dat <- split(as.matrix(dat), row(dat))
  lapply(dat, dat_to_mat)
}

# Calculate mean matrix by element - not a "real" matrix though: would median by-matrix (rather than by-element) 
# be better?

mat_mean <- function(mats) {
  n <- length(mats)
  mats <- array(unlist(mats), c(8, 8, n))
  rowMeans(mats, dims = 2)
}

# Same as above but for SE: standard error of the matrices across seeds

mat_se <- function(mats) {
  m <- length(mats)
  mats <- array(unlist(mats), c(8, 8, m))
  apply(mats, 1:2, std.error)
}


# Ideal way to use these functions: subset data by model index (or groups of them) and generation, 
# construct G matrices for each replicate via mat_gen
# then construct mean G matrix and CIs around them of each

# E.g. Model Index 2, gen 150,000

G_null_mat_mi2 <- mat_mean(mat_gen(d_null_mat[1:3,]))
G_null_mat_mi6 <- mat_mean(mat_gen(d_null_mat[4:6,]))
G_null_mat_mi8 <- mat_mean(mat_gen(d_null_mat[7:9,]))

G_null_mat <- array(data = c(G_null_mat_mi2, G_null_mat_mi6, G_null_mat_mi8), c(8, 8, 3))

# Eigentensor analysis: Using evolqg EigenTensorDecomposition()

library(evolqg)

EigenTensorDecomposition(G_null_mat)


# RandomSkewers analysis: also using evolqg, RandomSkewers()

G_null_mat <- list(G_null_mat_mi2, G_null_mat_mi6, G_null_mat_mi8)

RandomSkewers(G_null_mat)



