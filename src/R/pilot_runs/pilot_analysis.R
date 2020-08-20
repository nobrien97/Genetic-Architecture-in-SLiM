# Pilot data analysis: Take G matrices, eigenanalysis, tensors? 

setwd("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs")

# Basic standard error function for use with data manipulation

std.error <- function(x){
  n <- length(x)
  sd <- sd(x)
  sd/sqrt(n)
}


d_null_32bit <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/nimrod_tests/Output/32bit_seed/out_8T_null_means.csv", header = F)

# Name columns

names(d_null_32bit)[1:7]<- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu", "nloci")


# pleiocov terms                
names(d_null_32bit)[8:35] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))


names(d_null_32bit)[36] <- "Ne"

names(d_null_32bit)[37:44] <- paste0("mean", 0:7)

names(d_null_32bit)[45:52] <- paste0("var", 0:7)

names(d_null_32bit)[53:80] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null_32bit)[81:108] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null_32bit)[109:118] <- paste0("H_chr", 0:9)

d_null_32bit$seed <- as.factor(d_null_32bit$seed)



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

d_null$seed <- as.factor(d_null$seed)

# Organise data for G matrix construction

library(dplyr)

d_sbstfinal_null <- as_tibble(d_null[d_null$gen == 150000,])

d_null_mat <- select(as_tibble(d_null), gen, seed, modelindex, paste0("var", 0:7), c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7)))




# Data coercion to matrix function: expects a single row of 39 variables: gen, seed, modelindex, 8 variances, 28 covariances

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
  mode(M) = "numeric"
  M
}

# Split data frame up into single lines, do dat_to_mat on each, store as an array
# Split by seed first so tree is list -> gen -> seed -> model
# Easier to compare models this way by seed


mat_gen <- function(dat) {
  dat <- dplyr::arrange(dat, gen, modelindex, seed)
  dat <- group_split(dat, gen)
  dat <- lapply(dat, function(x) { dplyr::group_split(x, seed)})
  dat <- lapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y))
   })
  })
  dat <- lapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mat(z) 
      })
    }) 
  })
  dat
}

# Split by generation tests
dattest <- arrange(d_null_mat, gen, modelindex, seed)
dattest <- group_split(dattest, gen)
dattest <- lapply(dattest, function(x) { group_split(x, modelindex)})

dattest <- lapply(dattest, function(x) { lapply(x, function(y) {
  split(as.matrix(y), row(y))
  })
})
dattest <- lapply(dattest, function(x) { 
  lapply(x, function(y) { 
    lapply(y, function(z) { 
      dat_to_mat(z) 
      })
    }) 
  })


matgentest <- mat_gen(d_null_mat)

# Calculate mean matrix by element - not a "real" matrix though: would median by-matrix (rather than by-element) 
# be better?

mat_mean <- function(mats) {
  n <- length(mats)
  mats <- array(as.numeric(unlist(mats)), c(8, 8, n))
  rowMeans(mats, dims = 2)
}

# Same as above but for SE: standard error of the matrices across seeds

mat_se <- function(mats) {
  m <- length(mats)
  mats <- array(unlist(mats), c(8, 8, m))
  apply(mats, 1:2, std.error)
}

# Function to combine these so we put in a full data set, we split it into seeds, calculate a mean for each model
# and put them all in a list

matmean_construct <- function(dat) {
  matlist <- mat_gen(dat) # mat_gen splits each model into its seed rows and generation, stores them as a list
  lapply(matlist, function(x) {
    lapply(x, mat_mean) # apply mat_mean to each model
  })
}

matstruc_test <- matmean_construct(d_null_mat)

# Function to do as above, but for an se instead of mean

matse_construct <- function(dat) {
  matlist <- mat_gen(dat) # mat_gen splits each model into its seed rows, stores them as a list
  lapply(matlist, function(x) {
    lapply(x, mat_se)
    }) # apply mat_se to each model
}

## Multicore versions of the above functions
# Uses parallel, dplyr

MCmat_gen <- function(dat, cores) {
  dat <- dplyr::arrange(dat, gen, modelindex, seed)
  dat <- dplyr::group_split(dat, gen)
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed)}, mc.cores = cores)
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y))
  })
  }, mc.cores = cores)
  
  dat <- parallel::mclapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mat(z) 
      })
    }) 
  }, mc.cores = cores)
  dat
}

library(evolqg)
MCG_ET <- function(G, cores, nmats) {
  Gmax <- parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      mats <- simplify2array(y)
      evolqg::EigenTensorDecomposition(mats, return.projection = F)$matrices[1:nmats]
    })
  }, mc.cores = cores)
  Es <- parallel::mclapply(Gmax, function(x) {
    lapply(x, function(y) {
      eigen(y, symmetric = T, only.values = T)
    })
  }, mc.cores = cores)
  Es
}


is.positive.definite(testGs[[1]][[1]][,,1])
# first is subsetting by generation, then seed, then last is getting each model's matrix as element of array

G_ET <- function(G, nmats) {
  Gmax <- lapply(G, function(x) {
    lapply(x, function(y) {
      mats <- simplify2array(y)
      mats <- mats[,,matrixcalc::is.positive.definite(y) == T]
       evolqg::EigenTensorDecomposition(mats, return.projection = F)$matrices[1:nmats]
  #   })
  # })
  # Es <- lapply(Gmax, function(x) {
  #   lapply(x, function(y) {
  #     eigen(y, symmetric = T, only.values = T)
  #   })
  # })
  # Es
    }) })
  Gmax
}



testGs <- lapply(testmatrices, function(x) {
  lapply(x, function(y) {
      lapply(x, `[[`, 1)
    })
})
  
  n <- length(x)
  mats <- array(unlist(mats), c(8, 8, n))
  evolqg::EigenTensorDecomposition(mats, return.projection = F)$matrices[,,1:2]
})


library(evolqg)
mats_test <- array(unlist(testmatrices[[1]][1]), c(8,8,3))
ETOutput <- EigenTensorDecomposition(mats_test, return.projection = F)$matrices[,,1:2]
EigOutput <- eigen(ETOutput$matrices[2], symmetric = T, only.values = T)


# E.g. with d_null_mat

G_null_mean <- simplify2array(matstruc_test)

G_null_se <- as.array(matse_construct(d_null_mat))




# Eigentensor analysis: Using evolqg EigenTensorDecomposition()

library(evolqg)

Gmean_ET <- function(G) {
  Gmax <- lapply(G, function(x) {
      evolqg::EigenTensorDecomposition(x, return.projection = F)$matrices[,,1:2]
  })
    Es <- lapply(Gmax, function(x) {
      eigen(x, symmetric = T, only.values = T)
    })
    Es
}


Gmax_test <- lapply(matstruc_test, function(x) {
  evolqg::EigenTensorDecomposition(x, return.projection = F)$matrices[,,1:2]
})




ET_Decomp_Gnull <- EigenTensorDecomposition(G_null_mean)

E_Decomp_ETnull <- eigen(ET_Decomp_Gnull$matrices[,,1])

#> sum(E_Decomp_ETnull$values[1:2]) / sum(E_Decomp_ETnull$values)
#[1] 0.3819765
# Linear combinations of traits 1 and 2 of eigentensor 1 explain 38.2% of the divergence between populations (?) 



# RandomSkewers analysis: also using evolqg, RandomSkewers()

RandomSkewers(as.list(G_null_mat))

# rPCA
library(vcvComp)
eigen(mat.sq.dist(G_null_mean, dist. = "Euclidean"))


# VFA
library(Matrix)
#library(rARPACK)
library(vsp)

vsp()


# Recom model, same stuff



# recom data set

d_recom <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/pilot_runs/Output/out_8T_recom_means.csv", header = F)

# Name columns

names(d_recom)[1:7]<- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu", "nloci")


# pleiocov terms                
names(d_recom)[8:35] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))


names(d_recom)[36] <- "Ne"

names(d_recom)[37:44] <- paste0("mean", 0:7)

names(d_recom)[45:52] <- paste0("var", 0:7)

names(d_recom)[53:80] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_recom)[81:108] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_recom)[109:118] <- paste0("H_chr", 0:9)

d_recom$seed <- as.factor(d_recom$seed)

# Organise data for G matrix construction

library(dplyr)

d_sbstfinal_recom <- as_tibble(d_recom[d_recom$gen == 150000,])

d_recom_mat <- select(d_sbstfinal_recom, gen, seed, modelindex, paste0("var", 0:7), c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7)))


