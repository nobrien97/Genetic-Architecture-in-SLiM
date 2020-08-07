##############################################################################################################
#  Simple test script for including two traits, with arguments coming from latin square, for use on cluster
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Get environment variables

TMPDIR <- Sys.getenv('TMPDIR')
USER <- Sys.getenv('USER')


# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


# Seed generation - 100 seeds from the total range of possible values in a 64 bit signed integer

rsample <- as.character(runif(100, 1, (2^62 - 1)))


# Load LHS samples

lscombos <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/lscombos.csv"), header = T)

# - run on all available cores, treat them as nodes

cl <- makeCluster(future::availableWorkers(methods = "PBS"))
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to LHS in command line

foreach(i=rsample) %:%
  foreach(j=1:nrow(lscombos)) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -d seed=%i -d rregion=%s -d Ne=%i -d QTL_cov=%f -d pleiorate=%f -d delmu=%f -d rwide=%s -d modelindex=%i /home/$USER/SLiM/Scripts/neutral2T10GRjune_restruct2_opt.slim", 
                               as.character(i), as.character(lscombos$rregions[j]), as.integer(round(lscombos$Ne[j])), lscombos$pleiocov[j], lscombos$pleiorate[j], lscombos$delmu[j], as.character(lscombos$rwide[j]), j, intern=T))
  }
stopCluster(cl)
