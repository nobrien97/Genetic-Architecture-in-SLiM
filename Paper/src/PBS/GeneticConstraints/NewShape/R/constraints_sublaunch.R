##############################################################################################################
#  Test of the genetic constraints system: 24 replicates, stabilising selection, 4 genes, 3 QTLs
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Environment variables

USER <- Sys.getenv('USER')

# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


seeds <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Tools/SeedGen/seeds_sub.csv"), header = T)


cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM

foreach(i=seeds$Seed) %dopar% {
	# Use string manipulation functions to configure the command line args, feeding from a data frame of seeds
	# then run SLiM with system(),
    	slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s /home/$USER/SLiM/Scripts/Tests/GeneticConstraints/NewShape/slim/polygen_maint.slim", as.character(i), intern=T))
  }
stopCluster(cl)

