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

# mat_gen generates G matrices in nested list format - need to convert the last bit into arrays (so for each seed we have an array of matrices, one for each model index)
# Requires dplyr

library(dplyr)

mat_gen <- function(dat) {
  dat <- dplyr::arrange(dat, gen, modelindex, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen))
  dat <- lapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))})
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

test_datmat <- arrange(d_null_mat, gen, modelindex, seed)
test_datmat <- group_split(test_datmat, gen) %>% setNames(unique(test_datmat$gen))
test_datmat <- lapply(test_datmat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))})


# Test it on the data to make sure it works
matgentest <- mat_gen(d_null_mat)



## Multicore versions of the above functions
# Uses parallel, dplyr

MCmat_gen <- function(dat, cores) {
  dat <- dplyr::arrange(dat, gen, modelindex, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(y$modelindex)) # split the dataframe of seed values by their row (models), treat as a matrix instead of dataframe for dat_to_mat
  })
  }, mc.cores = cores)
  
  dat <- parallel::mclapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mat(z) # Convert each line (model/seed/generation combination) into a G matrix
      })
    }) 
  }, mc.cores = cores)
 
  dat
}

matgentest <- MCmat_gen(d_null_mat, 4)



# Ensure matrices are positive definite, get rid of any that aren't

MCG_PD <- function(G, cores) {
  G <- parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
        if (!matrixcalc::is.positive.definite(z))
          z <- NULL # Set the values that aren't positive definite as NULL to tag for removal later
        else
          z <- z # Required or the rest of the values are empty
      })
    })
  }, mc.cores = cores)
  G
}

matgentest_pdt <- MCG_PD(matgentest, 4)


# Eigentensor analysis - Take G matrices from mat_gen, make sure they are positive definite, decompose, ET analysis
# Using evolqg EigenTensorDecomposition()
# G is a list of lists of lists: gen -> seed -> models, each of these has a G matrix
# cores is number of cores to use concurrently when calculating eigentensors and eigenvalues
# nmats for the number of eigentensors to keep and do an eigenanalysis on
library(evolqg)
MCG_ET <- function(G, cores, nmats) {
  Gmax <- parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      eigen(evolqg::EigenTensorDecomposition(simplify2array(compact(y)), # Compact removes NULL values
                                             return.projection = F)$matrices[,,1:nmats], 
                                            symmetric = T, only.values = F)
    })
  }, mc.cores = cores)
  
 Gmax
}

testGs <- MCG_ET(matgentest_pdt[1], 4, 1)     # Just do it for one time point for trial



# RandomSkewers analysis: also using evolqg, RandomSkewers()


# Function for calculating random skewers between models

MCskew <- function(G, cores) {
  G <- parallel::mclapply(G, function(x) {
    lapply(x, function(y) { # RandomSkewers() act on list of matrices, cov.x
        evolqg::RandomSkewers(compact(y))
    })
  }, mc.cores = cores)
  G
}

testGskew <- MCskew(matgentest_pdt, 4)

# rPCA
library(vcvComp)
eigen(mat.sq.dist(G_null_mean, dist. = "Riemannian"))

# Multicore function to do Principal Coordinates Analysis on our data: adapted from https://cran.r-project.org/web/packages/vcvComp/vignettes/vcvComp-worked-example.html

MC_PCOA <- function(G, cores) {
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      sqdist <- vcvComp::mat.sq.dist(simplify2array(compact(y)), dist. = "Riemannian")
      ord <- vcvComp::pr.coord(sqdist)
      ord
    })
  }, mc.cores = cores)
}

# Test the rPCA function
G_rPCA <- MC_PCOA(matgentest_pdt, 4)

# Plot principal coordinates from PCOA: https://cran.r-project.org/web/packages/vcvComp/vignettes/vcvComp-worked-example.html

# Mean data across seeds for each model - need to flip list around to model first then seed?

MC_PCOA_plot <- function(G, cores) {
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      sqdist <- vcvComp::mat.sq.dist(simplify2array(compact(y)), dist. = "Riemannian")
      ord <- vcvComp::pr.coord(sqdist)$PCoords[,1:3]
      lapply(seq_len(nrow(ord)), function(z) {
        ord[z,]
      })
    })
  }, mc.cores = cores)
  
}



G_PCOA_plot <- MC_PCOA_plot(matgentest_pdt, 4)


d_null_arr <- arrange(d_null_mat %>% distinct(seed, gen, modelindex, .keep_all = T), gen, modelindex, seed)

test_df <- data.frame( # matgentest_pdt[[1]] is seed, compact([[1]][[1]]) is modelindex which is PD
  gen = rep(as.integer(names(matgentest)), each = length(matgentest_pdt[[1]])*length(compact(matgentest_pdt[[1]][[1]]))),
  seed = rep(names(matgentest[[1]]), each = length(compact(matgentest_pdt[[1]][[1]]))),
  modelindex = names(compact(matgentest_pdt[[1]][[1]])),
  pc <- paste(unlist(G_PCOA_plot), sep=",")
)

library(tidyverse)
G_PCOA_plot %>%
  map_df(transpose) %>%
  mutate_at(vars(c('PCo1', 'PCo2', 'PCo3')), funs(unlist))

# Visualization of PCo1, PCo2 and PCo3 (fig. 2)
coul.pop <- "blue" # colors
pco3d <- c(1, 2, 3)  # dimensions
xyzlab <- c(paste("PCoord", pco3d[1]), 
            paste("PCoord", pco3d[2]), 
            paste("PCoord", pco3d[3]))
s3d <- scatterplot3d:::scatterplot3d(G_rPCA[["150000"]][[1]]$PCoords[, pco3d[1:3]],
                                     xlab = xyzlab[1], ylab = xyzlab[2], zlab = xyzlab[3],
                                     color = coul.pop, pch = 19, angle = 55,
                                     type = "h", lty.hplot = 3, 
                                     cex.symbols = 1, cex.axis = 0.8)
s3d.coords <- s3d$xyz.convert(G_rPCA[["150000"]][[1]]$PCoords[, pco3d[1:3]])
text(s3d.coords$x, s3d.coords$y, 
     labels = row.names(G_rPCA[["150000"]][[1]]$PCoords), 
     pos = 4, cex = 0.7, col = coul.pop)






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



###################################
#   Deprecated - code graveyard   #
###################################

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

# Split by generation tests: used for making the function
dattest <- dplyr::arrange(d_null_mat, gen, modelindex, seed) # rearrange data so it is ordered when we list it
dattest <- dplyr::group_split(dattest, gen) # Split dataframe according to group
dattest <- lapply(dattest, function(x) { dplyr::group_split(x, modelindex)}) # split within that according to modelindex

dattest <- lapply(dattest, function(x) { lapply(x, function(y) {
  split(as.matrix(y), row(y)) # create a matrix from everything nested within that (individual seed/modelindex G matrices), splitting each row in the data frame into a separate one
})
})
dattest <- lapply(dattest, function(x) { 
  lapply(x, function(y) { 
    lapply(y, function(z) { 
      dat_to_mat(z) # convert each line within the list to a G matrix
    })
  }) 
})


n <- length(x)
mats <- array(unlist(mats), c(8, 8, n))
evolqg::EigenTensorDecomposition(mats, return.projection = F)$matrices[,,1:2]
})


# Combine the models into an array
testGs <- lapply(testmatrices, function(x) {
  lapply(x, function(y) {
    lapply(x, `[[`, 1) # Extracts values - ??? 
  })
})



# To return a single matrix from the list - e.g. at generation 50000, first seed, first model:
matgentest[["50000"]][[1]][[1]]
# or
matgentest[[1]][[1]][[1]]





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



is.positive.definite(testGs[[1]][[1]][[1]])
# first is subsetting by generation, then seed, then last is getting each model's matrix as element of array



library(evolqg)
mats_test <- array(unlist(testmatrices[[1]][1]), c(8,8,3))
ETOutput <- EigenTensorDecomposition(mats_test, return.projection = F)$matrices[,,1:2]
EigOutput <- eigen(ETOutput[,,1], symmetric = T, only.values = T)


# E.g. with d_null_mat

G_null_mean <- simplify2array(matstruc_test)

G_null_se <- as.array(matse_construct(d_null_mat))


# Function for calculating ET of mean matrix across seeds

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

# Single core test of eigentensors for error output/debug

G_ET <- function(G, nmats) {
  G <- lapply(G, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) { # make sure we only use the mats that are positive definite
        if (!matrixcalc::is.positive.definite(z)) 
        {
          z <- NULL
        }
        else 
          z <- z
      })
    })
  })
  GET <- lapply(G, function(x) {
    lapply(x, function(y) {
      y <- compact(simplify2array(y))
      evolqg::EigenTensorDecomposition(y, return.projection = F)$matrices[,,1:nmats]
    })
  })
  Es <- lapply(GET, function(x) {
    lapply(function(y) {
      eigen(y, symmetric = T, only.values = T)
    })
  })
  Es
}

# Test for diagnostics
testGs <- G_ET(matgentest, 1)


# Random Skewers test
RandomSkewers(as.list(G_null_mat))

# rPCA test

dists <- lapply(matgentest_pdt, function(x) {
  lapply(x, function(y) {
    G <- simplify2array(compact(y))
    vcvComp::mat.sq.dist(G, dist. = "Riemannian")
  })
})


test_PCOA <- lapply(matgentest_pdt[1:2], function(x) {
  lapply(x, function(y) {
    sqdist <- vcvComp::mat.sq.dist(simplify2array(compact(y)), dist. = "Riemannian")
    ord <- vcvComp::pr.coord(sqdist)$PCoords[,1:3]
    ord
    lapply(seq_len(nrow(ord)), function(z) ord[z,])
  })
})

sqdist <- vcvComp::mat.sq.dist(simplify2array(compact(matgentest_pdt[[1]][[1]])), dist. = "Riemannian")
ord <- vcvComp::pr.coord(sqdist)$PCoords[,1:3]

ord[1,]