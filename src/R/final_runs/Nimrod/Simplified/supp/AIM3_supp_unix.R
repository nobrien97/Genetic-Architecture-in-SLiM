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

lscombos_sel <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Final/AIM3/Inputs/lscombos_sel.csv"), header = T)

combo <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Final/AIM3/Inputs/AIM3_supp_combos.csv"), header = T)

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to combo - only need to run 24 models


# 5 Parameters for Aim 1, 6 for Aim 3 (tau):
# rwide, pleio_cov, pleiorate, delmu, locisigma, tau


foreach(i=combo$X) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
	j <- combo[i,]$model
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d locisigma=%f -d tau=%f -d modelindex=%i /home/$USER/SLiM/Scripts/Nimrod/Final/AIM3/Inputs/stabsel_recom_8T100L.slim", 
                               as.character(combo[i,]$seed), as.character(lscombos_sel[j,]$rwide), lscombos_sel[j,]$pleiocov, lscombos_sel[j,]$pleiorate, lscombos_sel[j,]$delmu, lscombos_sel[j,]$locisigma, ls_combos_sel[j,]$tau, j, intern=T))
  }

stopCluster(cl)
