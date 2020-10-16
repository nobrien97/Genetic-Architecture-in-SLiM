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
  means <- dat[39:46]
  means
}

# Function to rearrange data into list sorted by variables
mean_gen <- function(dat) {
  dat <- dplyr::arrange(dat, gen, modelindex, seed) # should be modelindex instead of tau in final function
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen))
  dat <- lapply(dat, function(x) { dplyr::group_split(x, seed)  %>% setNames(unique(x$seed))}) # Don't need names of seeds
  dat <- lapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(dat$modelindex))
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
opt_gen <- function(opt, v) {
  opt <- dplyr::arrange(opt, v, seed) # should be modelindex instead of tau in final function
  opt <- dplyr::group_split(opt, seed)  %>% setNames(sort(unique(opt$seed)))
  opt <- lapply(opt, function (x) {
    x <- x[,3:10]
    split(as.matrix(x), row(x)) %>% setNames(sort(unique(v)))
  })
  opt
}


# Euclidean distance: needs to be for each time point, model and seed combo

# Function needs to go into list, lapply level 1 is generation, level 2 is seed, level 3 is model/tau
# So go to level 3, calculate distance between that and the optimum (where optimum value is from another data frame)

MCeuc_dist <- function(dat, opt, cores) {
  require(parallel)
  dat <- mean_gen(dat)
  opt <- opt_gen(opt, opt$modelindex)
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
  combos <- sample(1:length(G[[1]][[1]]), 2*n, replace = T) # pairwise combinations, in a matrix  
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


MC_relW_PW_combn <- function(G, n=100, cores) {
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

# Relative PCA 3: comparing null against sel models
# Need to pick similar values for our tested interaction (e.g. bins of low recom + high delmu etc.), with different levels of selection
# So testing difference from NULL with logGV
# Will require different matrix structure with multiple predictors

# 1) Add an interaction term to use as a subsetter
# 2) Randomly choose a row which fits each unique interaction term for each tau category
# 3) Convert to matrix and run three rPCAs, one for each selection term (null vs low, null vs med, null vs high)
# 4) Store output in nested list
# 5) Repeat for number of replicates
# 6) Name list properly

rPCA <- function(dat, var1, var2, n) {
  require(dplyr)
  require(vcvComp)
  
  dat$interact <- interaction(var1, var2) # Interaction between the two

  # For each unique interaction, randomly select n rows with that interaction from null, low, medium, and high selection
  # Then dat_to_mat those rows, and do three rPCAs: null, low, medium, high, and store that output in a list
  # so logGV_low <- rPCA_null.low$logGV or something
  i <- 1
  rel_ls <- vector("list", n)
  while (i < (n+1)) { # While we're still running replicates, do the for loop (eek - but seems to run alright, not too slow)
  for (inter in unique(dat$interact)) { 
    sampled_rows <- dat[1:4,] # Create new data frame to fill in with the same columns as the original data, we overwrite the original rows
    for (sel in unique(dat$tau.cat)) {
      sampled_rows[match(sel, unique(dat$tau.cat)),] <- slice_sample(subset(dat, interact == inter & tau.cat == sel), 1) # Randomly sample from the proper subset of data (certain interaction and selection strength)
    }
    sample_null <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Null",] )
    sample_high <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "High",] )
    sample_med <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Medium",] )
    sample_low <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Low",] )
    
    rel_low <- vcvComp::relative.eigen(sample_null, sample_low)
    rel_med <- vcvComp::relative.eigen(sample_null, sample_med)
    rel_high <- vcvComp::relative.eigen(sample_null, sample_high)
    
    # Store output in list
    rel_ls[[i]][[match(inter, unique(dat$interact))]] <- list(Low = rel_low,
                                                              Med = rel_med,
                                                              High = rel_high)
    }
    # At the end of the loop (before we repeat), name the element of the list we just did (i.e. the replicate)
    names(rel_ls[[i]]) <- unique(dat$interact)
    i <- i+1
  }
  rel_ls
}





# This function is slightly different to dat_to_mat() because the d_raw_mat used has more columns for each of the parameters, 
# so indexing vars and covs is different
dat_to_mat_rPCA <- function(dat) {
  dat <- as.vector(t(dat))
  vars <- dat[15:22]
  covs <- dat[23:50]
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

MCOrg_rPCA <- function(relPCs, cores) {
  require(tidyverse)
  parallel::mclapply(relPCs, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
        list(relGmax.val = z$relValues[1],
             relG2.val = z$relValues[2],
             relGV = z$relGV,
             logGV = z$logGV,
             distCov = z$distCov,
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


# Eigentensors: Includes null, low, medium, high, characterises variability between all of them
# Need to pick similar values for our tested interaction (e.g. bins of low recom + high delmu etc.), with different levels of selection
# So testing difference from NULL with first eigentensor's eigenvalues
# Will require different matrix structure with multiple predictors

# 1) Add an interaction term to use as a subsetter
# 2) Randomly choose a row which fits each unique interaction term for each tau category
# 3) Convert to matrix and run three eigentensor decompositions, one for each selection term (null vs low, null vs med, null vs high)
# 4) Run regular eigenanlysis on first eigentensor for each term
# 5) Store output in nested list
# 6) Repeat for number of replicates
# 7) Name list properly

eigentensor_G <- function(dat, var1, var2, n, cores) {
  require(plyr)
  require(dplyr)
  require(parallel)
  require(evolqg)
  
  dat$interact <- interaction(var1, var2) # Interaction between the two
  
  # For each unique interaction, randomly select n rows with that interaction from null, low, medium, and high selection
  # Then dat_to_mat those rows, and do ET decomp between the models: null, low, medium, high, and store that output in a list
  ET_ls <- vector("list", length(unique(dat$interact)))
  ET <-  parallel::mclapply((1:n), function(iter) {
    lapply(unique(dat$interact), function(x) {
      sampled_rows <- dat[1:4,] # Create new data frame to fill in with the same columns as the original data, we overwrite the original rows
      lapply(unique(dat$tau.cat), function(y, sampled_rows) {
        sampled_rows[match(y, unique(dat$tau.cat)),] <<- slice_sample(subset(dat, interact == x & tau.cat == y), 1) # Randomly sample from the proper subset of data (certain interaction and selection strength)
      }, sampled_rows = sampled_rows)
      
      sample_null <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Null",] )
      sample_high <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "High",] )
      sample_med <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Medium",] )
      sample_low <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Low",] )
      
      ETs <- evolqg::EigenTensorDecomposition(simplify2array(list(sample_null, sample_low, sample_med, sample_high)))
      
      eig_ETs <- eigen(ETs$matrices[,,1]) # Eigenanalysis on the first eigentensor
      # Store output in list
      ET_ls[[match(x, unique(dat$interact))]] <<- list(ET = ETs,
                                                       ET1_eig = eig_ETs)
    })
    # At the end of the loop (before we repeat), name the element of the list we just did (i.e. the replicate)
    names(ET_ls) <<- unique(dat$interact)
    ET_ls
    }, mc.cores = cores)
  
  ET
  
}
  

# Organise into list with 1 value each list element, so unlist works properly: separate relative eigenvalues/vectors 
# into columns cuts down to only Gmax and G2
MCOrg_ETeig <- function(G, cores) {
  require(tidyverse)
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      list(Gmax.val = y[[2]]$values[1],
           G2.val = y[[2]]$values[2],
           G3.val = y[[2]]$values[3], 
           G4.val = y[[2]]$values[4],
           G5.val = y[[2]]$values[5],
           G6.val = y[[2]]$values[6],
           G7.val = y[[2]]$values[7],
           G8.val = y[[2]]$values[8],
           Gtot.val = sum(y[[2]]$values),
           Gmaxprop.val = y[[2]]$values[1]/(sum(y[[2]]$values)),
           Gmax.vec1 = y[[2]]$vectors[1,1],
           Gmax.vec2 = y[[2]]$vectors[2,1],
           Gmax.vec3 = y[[2]]$vectors[3,1],
           Gmax.vec4 = y[[2]]$vectors[4,1],
           Gmax.vec5 = y[[2]]$vectors[5,1],
           Gmax.vec6 = y[[2]]$vectors[6,1],
           Gmax.vec7 = y[[2]]$vectors[7,1],
           Gmax.vec8 = y[[2]]$vectors[8,1],
           G2.vec1 = y[[2]]$vectors[1,2],
           G2.vec2 = y[[2]]$vectors[2,2],
           G2.vec3 = y[[2]]$vectors[3,2],
           G2.vec4 = y[[2]]$vectors[4,2],
           G2.vec5 = y[[2]]$vectors[5,2],
           G2.vec6 = y[[2]]$vectors[6,2],
           G2.vec7 = y[[2]]$vectors[7,2],
           G2.vec8 = y[[2]]$vectors[8,2],
           G3.vec1 = y[[2]]$vectors[1,3],
           G3.vec2 = y[[2]]$vectors[2,3],
           G3.vec3 = y[[2]]$vectors[3,3],
           G3.vec4 = y[[2]]$vectors[4,3],
           G3.vec5 = y[[2]]$vectors[5,3],
           G3.vec6 = y[[2]]$vectors[6,3],
           G3.vec7 = y[[2]]$vectors[7,3],
           G3.vec8 = y[[2]]$vectors[8,3],
           G4.vec1 = y[[2]]$vectors[1,4],
           G4.vec2 = y[[2]]$vectors[2,4],
           G4.vec3 = y[[2]]$vectors[3,4],
           G4.vec4 = y[[2]]$vectors[4,4],
           G4.vec5 = y[[2]]$vectors[5,4],
           G4.vec6 = y[[2]]$vectors[6,4],
           G4.vec7 = y[[2]]$vectors[7,4],
           G4.vec8 = y[[2]]$vectors[8,4],
           G5.vec1 = y[[2]]$vectors[1,5],
           G5.vec2 = y[[2]]$vectors[2,5],
           G5.vec3 = y[[2]]$vectors[3,5],
           G5.vec4 = y[[2]]$vectors[4,5],
           G5.vec5 = y[[2]]$vectors[5,5],
           G5.vec6 = y[[2]]$vectors[6,5],
           G5.vec7 = y[[2]]$vectors[7,5],
           G5.vec8 = y[[2]]$vectors[8,5],
           G6.vec1 = y[[2]]$vectors[1,6],
           G6.vec2 = y[[2]]$vectors[2,6],
           G6.vec3 = y[[2]]$vectors[3,6],
           G6.vec4 = y[[2]]$vectors[4,6],
           G6.vec5 = y[[2]]$vectors[5,6],
           G6.vec6 = y[[2]]$vectors[6,6],
           G6.vec7 = y[[2]]$vectors[7,6],
           G6.vec8 = y[[2]]$vectors[8,6],
           G7.vec1 = y[[2]]$vectors[1,7],
           G7.vec2 = y[[2]]$vectors[2,7],
           G7.vec3 = y[[2]]$vectors[3,7],
           G7.vec4 = y[[2]]$vectors[4,7],
           G7.vec5 = y[[2]]$vectors[5,7],
           G7.vec6 = y[[2]]$vectors[6,7],
           G7.vec7 = y[[2]]$vectors[7,7],
           G7.vec8 = y[[2]]$vectors[8,7],
           G8.vec1 = y[[2]]$vectors[1,8],
           G8.vec2 = y[[2]]$vectors[2,8],
           G8.vec3 = y[[2]]$vectors[3,8],
           G8.vec4 = y[[2]]$vectors[4,8],
           G8.vec5 = y[[2]]$vectors[5,8],
           G8.vec6 = y[[2]]$vectors[6,8],
           G8.vec7 = y[[2]]$vectors[7,8],
           G8.vec8 = y[[2]]$vectors[8,8]
       )
    })
  }, mc.cores = cores)
}

MCOrg_ET <- function(G, cores) {
  require(tidyverse)
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      VCV <- y[[1]]$matrices[,,1]
        list(ET1.val = y[[1]]$values[1],
             ET_tot.val = sum(y[[1]]$values),
             ET_prop.val = (y[[1]]$values[1])/(sum(y[[1]]$values)),
             ET1.proj_N = y[[1]]$projection[1,1],
             ET1.proj_L = y[[1]]$projection[2,1],
             ET1.proj_M = y[[1]]$projection[3,1],
             ET1.proj_H = y[[1]]$projection[4,1],
             ET1.mat_V0 = VCV[1,1],
             ET1.mat_V1 = VCV[2,2],
             ET1.mat_V2 = VCV[3,3],
             ET1.mat_V3 = VCV[4,4],
             ET1.mat_V4 = VCV[5,5],
             ET1.mat_V5 = VCV[6,6],
             ET1.mat_V6 = VCV[7,7],
             ET1.mat_V7 = VCV[8,8],
             ET1.mat_CV01 = VCV[1,2],
             ET1.mat_CV02 = VCV[1,3],
             ET1.mat_CV03 = VCV[1,4],
             ET1.mat_CV04 = VCV[1,5],
             ET1.mat_CV05 = VCV[1,6],
             ET1.mat_CV06 = VCV[1,7],
             ET1.mat_CV07 = VCV[1,8],
             ET1.mat_CV12 = VCV[2,3],
             ET1.mat_CV13 = VCV[2,4],
             ET1.mat_CV14 = VCV[2,5],
             ET1.mat_CV15 = VCV[2,6],
             ET1.mat_CV16 = VCV[2,7],
             ET1.mat_CV17 = VCV[2,8],
             ET1.mat_CV23 = VCV[3,4],
             ET1.mat_CV24 = VCV[3,5],
             ET1.mat_CV25 = VCV[3,6],
             ET1.mat_CV26 = VCV[3,7],
             ET1.mat_CV27 = VCV[3,8],
             ET1.mat_CV34 = VCV[4,5],
             ET1.mat_CV35 = VCV[4,6],
             ET1.mat_CV36 = VCV[4,7],
             ET1.mat_CV37 = VCV[4,8],
             ET1.mat_CV45 = VCV[5,6],
             ET1.mat_CV46 = VCV[5,7],
             ET1.mat_CV47 = VCV[5,8],
             ET1.mat_CV56 = VCV[6,7],
             ET1.mat_CV57 = VCV[6,8],
             ET1.mat_CV67 = VCV[7,8]
       )
   })
  }, mc.cores = cores)
}





ListToDF_ET <- function(Glist, responses) {
  require(rrapply)
  require(tidyr)
  rrapply(Glist, how = "melt") %>%
    pivot_wider(names_from = "L3") %>%
    unnest(responses)
}


# New Euclidean distance functions


# Euclidean distances: functions for transforming trait means and optima into appropriate lists for analysis

# Function to convert a single line from dataframe to just the mean values
dat_to_mean_new <- function(dat) {
  dat <- as.vector(t(dat))
  means <- dat[8:15]
  means
}


# Multithread get mean values from populations, store in list as usual
MCmean_gen_new <- function(dat, cores) {
  require(parallel)
  require(dplyr)
  dat <- dplyr::arrange(dat, gen, modelindex, seed)
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen)) # split data frame by generation
  dat <- parallel::mclapply(dat, function(x) { dplyr::group_split(x, seed) %>% setNames(unique(x$seed))}, mc.cores = cores) # split by seed
  dat <- parallel::mclapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(dat$modelindex)) # split the dataframe of seed values by their row (models), treat as a matrix instead of dataframe for dat_to_mat
  })
  }, mc.cores = cores)
  dat
  dat <- parallel::mclapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mean_new(z) # Convert each line (model/seed/generation combination) into a G matrix
      })
    }) 
  }, mc.cores = cores)
  
  dat
}

# Convert optimums into a list for easier comparison
MCopt_gen <- function(opt, cores) {
  require(parallel)
  opt <- dplyr::arrange(opt, modelindex, seed) # should be modelindex instead of tau in final function
  opt <- dplyr::group_split(opt, seed)  %>% setNames(sort(unique(opt$seed)))
  opt <- parallel::mclapply(opt, function (x) {
    x <- x[,3:10]
    split(as.matrix(x), row(x)) %>% setNames(sort(unique(v)))
  }, mc.cores = cores)
  opt
}

mean_gen_new <- function(dat) {
  dat <- dplyr::arrange(dat, gen, modelindex, seed) # should be modelindex instead of tau in final function
  dat <- dplyr::group_split(dat, gen) %>% setNames(unique(dat$gen))
  dat <- lapply(dat, function(x) { dplyr::group_split(x, seed)  %>% setNames(unique(x$seed))}) # Don't need names of seeds
  dat <- lapply(dat, function(x) { lapply(x, function(y) {
    split(as.matrix(y), row(y)) %>% setNames(unique(dat$modelindex))
  })
  })
  dat <- lapply(dat, function(x) { 
    lapply(x, function(y) { 
      lapply(y, function(z) { 
        dat_to_mean_new(z) 
      })
    }) 
  })
  dat
}


MCeuc_dist_new <- function(dat, opt, cores) {
  require(parallel)
  dat <- MCmean_gen_new(dat, cores)
  opt <- MCopt_gen(opt, cores)
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



