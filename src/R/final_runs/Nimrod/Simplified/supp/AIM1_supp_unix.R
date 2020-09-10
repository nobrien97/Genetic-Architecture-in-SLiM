##############################################################################################################
#  Final Null model script with arguments coming from latin square, for use on cluster
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


# Load LHS samples - We are making up for what wasn't run the first time, re-running with same seed/model combo

lscombos_null <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Final/AIM1/Inputs/lscombos_null.csv"), header = T)

combo <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Final/AIM1/Inputs/AIM1_supp_combos.csv"), header = T)

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to combo - only need to run 25 models


# 5 Parameters for Aim 1, 6 for Aim 3 (wsd):
# rwide, pleio_cov, pleiorate, delmu, locisigma


foreach(i=combo$modelindex) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d locisigma=%f -d modelindex=%i /home/$USER/SLiM/Scripts/Nimrod/Final/AIM1/Inputs/null8T100L.slim", 
                               as.character(combo[match(i, combo$modelindex),]$seed), as.character(lscombos_null[i,]$rwide), lscombos_null[i,]$pleiocov, lscombos_null[i,]$pleiorate, lscombos_null[i,]$delmu, lscombos_null[i,]$locisigma, i, intern=T))
  }

stopCluster(cl)
