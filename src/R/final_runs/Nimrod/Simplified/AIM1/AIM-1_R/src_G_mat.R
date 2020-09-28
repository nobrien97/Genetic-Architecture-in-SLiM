###############################################################
#           Functions to wrangle G matrices into list         #
###############################################################


# Basic standard error function for use with data manipulation

std.error <- function(x){
  n <- length(x)
  sd <- sd(x)
  sd/sqrt(n)
}

# Sampling test: pairwise euclidean distances between scaled hypercube points
# latin hypercube homogeneity test

LHC.homogen <- function(LHC, samples, cores) {
  require(parallel)
  require(scales)
  dat <- LHC[,-1] 
  # rescale combos to between 0 and 1
  for (col in colnames(dat)) {
    dat[col] <- scales::rescale(dat[,col])
  }
  dat
  combos <- combn(sample(1:length(dat[,1]), samples), 2) # pairwise combinations, in matrix
  cmb_seq <- seq(1, length(combos), by = 2) # Comparisons between elements x and x+1 in combos
  d_dist <- double(length = 2*length(combos))
  
  d_dist_list <- mclapply(cmb_seq, function(x) {
                    x1 <- c(t(dat[combos[x],]))
                    x2 <- c(t(dat[combos[x+1],]))
                    d_dist[x] <- dist(rbind(x1, x2))
                    d_dist
                  }, mc.cores = cores)
  d_dist_list <- unlist(d_dist_list)
  d_dist <- d_dist_list[which(d_dist_list > 0.0)]
  range(d_dist)
}


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

# 2 trait version

dat_to_mat_2T <- function(dat) {
  dat <- as.vector(t(dat))
  vars <- dat[4:5]
  covs <- dat[12]
  M <- matrix(data = c(
    vars[1], covs[1],
    covs[1], vars[2]), nrow = 2, byrow = T
  )
  mode(M) = "numeric"
  M
}


# mat_gen generates G matrices in nested list format - need to convert the last bit into arrays (so for each seed we have an array of matrices, one for each model index)
# Requires dplyr
# Uses parallel, dplyr

MCmat_gen <- function(dat, v, cores) {
  require(parallel)
  require(dplyr)
  dat <- dplyr::arrange(dat, gen, v, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(v)) # split the dataframe of seed values by their row (models), treat as a matrix instead of dataframe for dat_to_mat
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

# 2T version

MCmat2T_gen <- function(dat, v, cores) {
  require(parallel)
  require(dplyr)
  dat <- dplyr::arrange(dat, gen, v, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(v)) # split the dataframe of seed values by their row (models), treat as a matrix instead of dataframe for dat_to_mat
  })
  }, mc.cores = cores)
  
  dat <- parallel::mclapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mat_2T(z) # Convert each line (model/seed/generation combination) into a G matrix
      })
    }) 
  }, mc.cores = cores)
  
  dat
}


# Version in opposite order, model first then seed

MCmatMS_gen <- function(dat, v, modnames, cores) {
  require(parallel)
  require(dplyr)
  group <- enquo(v)
  dat <- dplyr::arrange(dat, gen, v, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, !!group) %>% setNames(modnames)}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(y$seed) # split the dataframe of models by their row (seed), treat as a matrix instead of dataframe for dat_to_mat
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


# Ensure matrices are positive definite, set every one that isn't to NULL, to be removed as needed with compact()

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

# Eigentensor analysis - Take G matrices from mat_gen, make sure they are positive definite, decompose, ET analysis
# Using evolqg EigenTensorDecomposition()
# G is a list of lists of lists: gen -> seed -> models, each of these has a G matrix
# cores is number of cores to use concurrently when calculating eigentensors and eigenvalues
# nmats for the number of eigentensors to keep and do an eigenanalysis on
MCG_ET <- function(G, cores, nmats) {
  require(evolqg)
  Gmax <- parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      eigen(evolqg::EigenTensorDecomposition(simplify2array(compact(y)), # Compact removes NULL values
                                             return.projection = F)$matrices[,,1:nmats], 
            symmetric = T, only.values = F)
    })
  }, mc.cores = cores)
  
  Gmax
}


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

# Convert PCOA list to a data frame

PCOA_df <- function(PCOA) {
  require(dplyr)
  PCOA %>%
    map_df(transpose) %>%
    mutate_at(vars(c('PCo1', 'PCo2', 'PCo3')), funs(unlist))
}




# Euclidean distances: functions for transforming trait means and optima into appropriate lists for analysis

# Function to convert a single line from dataframe to just the mean values
dat_to_mean <- function(dat) {
  dat <- as.vector(t(dat))
  means <- dat[35:42]
  means
}

# Function to rearrange data into list sorted by variables
mean_gen <- function(dat) {
  dat <- dplyr::arrange(dat, gen, tau, seed) # should be modelindex instead of tau in final function
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen))
  dat <- lapply(dat, function(x) { dplyr::group_split(x, seed) }) # %>% setNames(unique(x$seed))}) Don't need names of seeds
  dat <- lapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y))
  })
  })
  dat <- lapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mean(z) 
      })
    }) 
  })
  dat
}


# Multithread get mean values from populations, store in list as usual
MCmean_gen <- function(dat, v, cores) {
  require(parallel)
  require(dplyr)
  dat <- dplyr::arrange(dat, gen, v, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(v)) # split the dataframe of seed values by their row (models), treat as a matrix instead of dataframe for dat_to_mat
  })
  }, mc.cores = cores)
  dat
  dat <- parallel::mclapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mean(z) # Convert each line (model/seed/generation combination) into a G matrix
      })
    }) 
  }, mc.cores = cores)
  
  dat
}

# Convert optimums into a list for easier comparison
opt_gen <- function(opt) {
  opt <- dplyr::arrange(opt, tau, seed) # should be modelindex instead of tau in final function
  opt <- dplyr::group_split(opt, seed) # %>% setNames(unique(opt$seed)) Don't need names of seeds
  opt <- lapply(opt, function (x) {
    x <- x[,3:10]
    split(as.matrix(x), row(x))
  })
  opt
}


# Euclidean distance: needs to be for each time point, model and seed combo

# Function needs to go into list, lapply level 1 is generation, level 2 is seed, level 3 is model/tau
# So go to level 3, calculate distance between that and the optimum (where optimum value is from another data frame)

MCeuc_dist <- function(dat, opt, cores) {
  require(parallel)
  dat <- mean_gen(dat)
  opt <- opt_gen(opt)
  dists <- parallel::mclapply(seq_along(dat), function(x) {
    lapply(seq_along(dat[[x]]), function(y) {
      lapply(seq_along(dat[[x]][[y]]), function(z) {
        opt_x <- as.numeric(opt[[y]][[z]]) # Get the index of the second and third levels of the list (seed, modelindex); for use in getting the right opt
        dist_x <- as.numeric(dat[[x]][[y]][[z]]) # Choose the right sampled vector of means (gen, seed, modelindex)
        dist(rbind(dist_x, opt_x)) # Euclidean distance calculation between given vector of means and associated optimum vector
      })
    })
  }, mc.cores = cores)
  dists
}


# Convert to data frame

dist_df <- function(dat, dists) {
  data.frame(
    gen = rep(unique(dat$gen), each = length(unique(dat$seed))*length(unique(dat$modelindex))),
    seed = rep(unique(dat$seed), each = length(unique(dat$modelindex))),
    modelindex = unique(dat$modelindex),
    distance = unlist(dists)
  )
}
 

## Regular old PCA on each model

MCPCA <- function(G, cores) {
  require(tidyverse)
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
        ev <- eigen(z, symmetric = T)
        as.list(ev)
      })
    })
  }, mc.cores = cores)
  
}

# Organise into list with 1 value each list element, so unlist works properly: separate eigenvalues/vectors into columns
# Makes S3 eigen into list, cuts down to only Gmax and G2
MCOrg_G <- function(G, cores) {
  require(tidyverse)
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
          list(Gmax.val = z$values[1],
               G2.val = z$values[2],
               Gmax.vec1 = z$vectors[1,1],
               Gmax.vec2 = z$vectors[2,1],
               Gmax.vec3 = z$vectors[3,1],
               Gmax.vec4 = z$vectors[4,1],
               Gmax.vec5 = z$vectors[5,1],
               Gmax.vec6 = z$vectors[6,1],
               Gmax.vec7 = z$vectors[7,1],
               Gmax.vec8 = z$vectors[8,1],
               G2.vec1 = z$vectors[1,2],
               G2.vec2 = z$vectors[2,2],
               G2.vec3 = z$vectors[3,2],
               G2.vec4 = z$vectors[4,2],
               G2.vec5 = z$vectors[5,2],
               G2.vec6 = z$vectors[6,2],
               G2.vec7 = z$vectors[7,2],
               G2.vec8 = z$vectors[8,2]
               )
      })
    })
  }, mc.cores = cores)
}


# Organise means

MCOrg_means <- function(means, cores) {
  require(tidyverse)
  parallel::mclapply(means, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
        list(T0_mean = as.numeric(z[1]),
             T1_mean = as.numeric(z[2]),
             T2_mean = as.numeric(z[3]),
             T3_mean = as.numeric(z[4]),
             T4_mean = as.numeric(z[5]),
             T5_mean = as.numeric(z[6]),
             T6_mean = as.numeric(z[7]),
             T7_mean = as.numeric(z[8])
        )
      })
    })
  }, mc.cores = cores)
}


# Nested List to Data frame using rrapply
# Thanks to answers at: https://stackoverflow.com/questions/63895533/converting-a-deeply-nested-list-to-a-dataframe

ListToDF <- function(Glist, responses, predictor) {
  require(rrapply)
  require(tidyr)
  rrapply(Glist, how = "melt") %>%
    pivot_wider(names_from = "L4") %>%
    unnest(responses) %>%
    rename(gen = L1, seed = L2, {{ predictor }} := L3)
}


# Fst calculation

MCFst <- function(dat, v, cores) {
  require(parallel)
  require(dplyr)
  dat <- dplyr::arrange(dat, gen, v, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), (0:(nrow(y) %/% 2))) %>% setNames(unique(v)) # split the dataframe of seed values into sets of 2 rows (delmu), treat as a matrix instead of dataframe for dat_to_mat
  })
  }, mc.cores = cores)
  
  dat <- parallel::mclapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        Fst_calc(z) # Convert each line (model/seed/generation combination) into a G matrix
      })
    }) 
  }, mc.cores = cores)
  
  dat
}


# Don't call this directly, called from within Fst()
Fst_calc <- function(dat) {
  H <- dat[107]
  Ft <- H[1,] + H[2,]
  Fs <- H
  Fst <- (Fs - Ft)/ Ft
  Fst
}


# within-group relative PCA: since so many combinations, just randomly sample a certain number

MC_relW.multi <- function(G, n=100, cores) {
  require(parallel)
  require(vcvComp)
  G2cmp <- sample(1:length(G[[1]][[1]]), n) # Sample some matrices to use for comparison
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      ycmp <- y[G2cmp]
      vcvComp::relGV.multi(simplify2array(ycmp))
    })
  }, mc.cores = cores)
}

# within-group relative PCA 2: this time with relative.eigen and eigen.test

MC_relW <- function(G, n=100, cores) {
  require(parallel)
  require(vcvComp)
  require(dplyr)
  if (n %% 2 > 0) # if n is odd, then add one so we can properly compare
    n <- n+1
  G2cmp <- sample(1:length(G[[1]][[1]]), n) # Sample some matrices to use for comparison
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      ycmp <- y[G2cmp]
      ycmpseq <- seq(1, length(ycmp), by = 2)
      lapply(ycmpseq, function(z) {
        mat1 <- matrix(unlist(ycmp[z]), nrow = 8)
        mat2 <- matrix(unlist(ycmp[z+1]), nrow = 8)
        vcvComp::relative.eigen(mat1, mat2)
      })
    })
  }, mc.cores = cores)
  
}

# Same as above, but with a different comparison method using combn to find pairwise combinations

MC_relW_PW <- function(G, n=100, cores) {
  require(parallel)
  require(vcvComp)
  require(dplyr)
  combos <- combn(sample(1:length(G[[1]][[1]]), n), 2) # pairwise combinations, in a matrix  
  lapply(G, function(x) {
    parallel::mclapply(x, function(y) {
      ycmpseq <- seq(1, length(combos), by = 2) # Comparisons between elements z and z+1 in combos
      lapply(ycmpseq, function(z) {
        name <- paste(names(y[combos[z]]), names(y[combos[z+1]]), sep = "_")
        mat1 <- matrix(unlist(y[combos[z]]), nrow = 8)
        mat2 <- matrix(unlist(y[combos[z+1]]), nrow = 8)
        eig <- vcvComp::relative.eigen(mat1, mat2)
        eig[["name"]] <- name
        eig
      })
    }, mc.cores = cores)
  })

}



# Organise into list with 1 value each list element, so unlist works properly: separate relative eigenvalues/vectors 
# into columns cuts down to only Gmax and G2
MCOrg_relG <- function(G, cores) {
  require(tidyverse)
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
         list(name = z$name,
             relGmax.val = z$relValues[1],
             relG2.val = z$relValues[2],
             logGV = z$logGV,
             relGmax.vec1 = z$relVectors[1,1],
             relGmax.vec2 = z$relVectors[2,1],
             relGmax.vec3 = z$relVectors[3,1],
             relGmax.vec4 = z$relVectors[4,1],
             relGmax.vec5 = z$relVectors[5,1],
             relGmax.vec6 = z$relVectors[6,1],
             relGmax.vec7 = z$relVectors[7,1],
             relGmax.vec8 = z$relVectors[8,1],
             relG2.vec1 = z$relVectors[1,2],
             relG2.vec2 = z$relVectors[2,2],
             relG2.vec3 = z$relVectors[3,2],
             relG2.vec4 = z$relVectors[4,2],
             relG2.vec5 = z$relVectors[5,2],
             relG2.vec6 = z$relVectors[6,2],
             relG2.vec7 = z$relVectors[7,2],
             relG2.vec8 = z$relVectors[8,2]
        )
      })
    })
  }, mc.cores = cores)
}
