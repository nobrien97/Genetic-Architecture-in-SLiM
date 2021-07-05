##############################################################################################################
#  Gadi single core speed test
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Environment variables

HOME <- Sys.getenv('HOME')

# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


seeds <- read.csv(paste0(HOME, "/Tools/SeedGenerator/seeds.csv"), header = T)


cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM

foreach(i=seeds$Seed) %dopar% {
	# Use string manipulation functions to configure the command line args, feeding from a data frame of seeds
	# then run SLiM with system(),
	slim_out <- system(sprintf("$HOME/SLiM/slim -s %s $HOME/SLiM/Tests/May2021_GadiSLiM_Speedtest/slim/polygen_maint.slim", as.character(i), intern=T))
  }
stopCluster(cl)

