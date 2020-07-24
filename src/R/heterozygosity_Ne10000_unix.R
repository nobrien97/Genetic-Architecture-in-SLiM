##############################################################################################################
#  Simple test script for including two traits, with arguments coming from latin square, for use on cluster
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.


TMPDIR <- Sys.getenv('TMPDIR')
USER <- Sys.getenv('USER')
PBSIND <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


# Seed generation - 100 seeds from the total range of possible values in a 64 bit signed integer

# rsample <- as.character(runif(100, 1, (2^62 - 1)))


# Load LHS samples - Split into 20 parts, so each node does 5 each (500 total runs/node)

rsample <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/seeds_het.csv"), header = F)
rsample <- as.character(as.vector(t(rsample)))

# - Multinode parallelism: Use the array index to tell us what range in lscombos to evaluate: from (index*5) -4: index*5

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to LHS in command line - j is from (PBSIND*5 )-4 to PBSIND*5


# Parameters:
# Ne, rregion, rwide, delmu, nloci (?), locisigma

Ne <- c(10000)
rwide <- c("0.0", "1.346e-5")
delmu <- c(0.0, 1.0)


foreach(i=rsample) %:%
  foreach(j=Ne) %:% 
    foreach(k=rwide) %:%
      foreach(m=delmu)  %dopar%   {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%s -d delmu=%f /home/$USER/SLiM/Scripts/heterozygosity_burnin_test.slim", 
                               as.character(i), j, k, m, intern=T))
  }
stopCluster(cl)
