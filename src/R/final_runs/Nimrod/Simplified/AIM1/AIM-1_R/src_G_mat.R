###############################################################
#           Functions to wrangle G matrices into list         #
###############################################################


# Basic standard error function for use with data manipulation

std.error <- function(x){
  n <- length(x)
  sd <- sd(x)
  sd/sqrt(n)
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
  means <- dat[36:43]
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
MCOrg <- function(G, cores) {
  require(tidyverse)
  parallel::mclapply(G, function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
          list(Gmax.val = z$values[1],
               G2.val = z$values[2],
               Gmax.vec = z$vectors[,1],
               G2.vec = z$vectors[,2])
      })
    })
  }, mc.cores = cores)
}


# https://stackoverflow.com/a/41882883
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}


