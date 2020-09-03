# Analysis ofdifferent Tau values and how they affect movement to the optimum over time

d_tau_means <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/Tau_test/Output/out_8T_stabsel_means.csv", header = F)

d_tau_opt <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/Tau_test/Output/out_8T_stabsel_opt.csv", header = F)

# Names

names(d_tau_means)[1:7] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu", "tau")

names(d_tau_means)[8:35] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_tau_means)[36:43] <- paste0("mean", 0:7)

names(d_tau_means)[44:51] <- paste0("var", 0:7)

names(d_tau_means)[52:79] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_tau_means)[80:107] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_tau_means)[108] <- "H"


names(d_tau_opt) <- c("seed", "tau", paste0("opt", 0:7))

# Join the frames together
library(dplyr)
d_master <- right_join(d_tau_opt[,2:10], d_tau_means, by = "tau")

#dists <- numeric(length = length(unique(modelindex))*length(unique(seed))* length(unique(gen)))
dists <- numeric(length = 261)
for (i in unique(d_tau_means$seed)) {
  for (j in unique(d_tau_means$tau)) {
    for (k in unique(d_tau_means$gen)) {
    dist_ex <- as.numeric(d_tau_means[d_tau_means$gen == unique(d_tau_means$gen)[k] & d_tau_means$tau == unique(d_tau_means$tau)[j] & d_tau_means$seed == unique(d_tau_means$seed)[i],][36:43])
    dists[] <- dist(rbind(dist_ex, as.numeric(d_tau_opt[d_tau_opt$seed == unique(d_tau_opt$seed)[i] & d_tau_opt$tau == unique(d_tau_opt$tau)[j],][3:10])))
    }
    
  }
}

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

d_tau_nodup <- d_tau_means %>% distinct(seed, gen, tau, .keep_all = T)
dat_test <- arrange(d_tau_means, gen, tau, seed)

# Test the functions

mean_list <- mean_gen(d_tau_nodup)

opt_list <- opt_gen(d_tau_opt)


# Euclidean distance: needs to be for each time point, model and seed combo

# Function needs to go into list, lapply level 1 is generation, level 2 is seed, level 3 is model/tau
# So go to level 3, calculate distance between that and the optimum (where optimum value is from another data frame)

euc_dist <- function(dat, opt) {
  dat <- mean_gen(dat)
  opt <- opt_gen(opt)
  dists <- lapply(seq_along(dat), function(x) {
    lapply(seq_along(dat[[x]]), function(y) {
      lapply(seq_along(dat[[x]][[y]]), function(z) {
        opt_x <- as.numeric(opt[[y]][[z]]) # Get the index of the second and third levels of the list (seed, modelindex); for use in getting the right opt
        dist_x <- as.numeric(dat[[x]][[y]][[z]]) # Choose the right sampled vector of means (gen, seed, modelindex)
        dist(rbind(dist_x, opt_x)) # Euclidean distance calculation between given vector of means and associated optimum vector
      })
    })
  })
  dists
}

euc_test <- euc_dist(d_tau_nodup, d_tau_opt)


# Convert to data frame for plotting: need to take outer list as gen, next level as seed, and next level as model





# Code Graveyard: Testing stuff

mean_list_lv1 <- mean_list[[1]][[1]]
mean_list_lv2 <- mean_list[[1]][[1]][[1]] # Checking subsets to decide how many levels we need to get to the right level in the data

match(mean_list_lv2, unlist(mean_list_lv1)) # Attempting to get the index of a lapply loop with match() - didn't work, and theoretically kind of slow so I abandoned that


opt <- arrange(d_tau_opt, tau, seed) # Double checking to see how arranging data would influence order of the nested list

# Example distance calculations for validation that this function ( dist() ) actually works as expected
dist_ex <- as.numeric(d_tau_means[d_tau_means$gen == 15000 & d_tau_means$tau == 10 & d_tau_means$seed == 167,][36:43])

dist_ex2 <- as.numeric(d_tau_means[d_tau_means$gen == 15000 & d_tau_means$tau == 1000 & d_tau_means$seed == 167,][36:43])

dist_ex3 <- as.numeric(d_tau_means[d_tau_means$gen == 15000 & d_tau_means$tau == 100 & d_tau_means$seed == 167,][36:43])


dist(rbind(dist_ex3, as.numeric(d_tau_opt[6,][3:10])))
