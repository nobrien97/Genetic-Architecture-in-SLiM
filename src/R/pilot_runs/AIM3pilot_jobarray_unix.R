##############################################################################################################
#  Pilot run stabilising selection model script with arguments coming from latin square, for use on cluster
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


# Load LHS samples - For this test run, we are doing 5 LHC samples and 3 seeds

lscombos_sel_pilot <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/lscombos_sel_pilot.csv"), header = T)

rsample <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/seeds_pilot.csv"), header = F)
rsample <- as.character(as.vector(t(rsample)))


cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to LHC samples in command line


# 8 Parameters for Aim 1 and 2, 9 for Aim 3 (wsd):
# Ne, rwide, pleio_cov, pleiorate, delmu, nloci, locisigma, delchr, wsd


foreach(i=rsample) %:%
  foreach(j=1:(nrow(lscombos_sel_pilot_pilot))) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d nloci=%i -d locisigma=%f -d delchr=%i -d wsd=%f -d modelindex=%i /home/$USER/SLiM/Scripts/stabsel_recom_8T100L.slim", 
                               as.character(i), as.integer(round(lscombos_sel_pilot$Ne[j])), as.character(lscombos_sel_pilot$rwide[j]), lscombos_sel_pilot$pleiocov[j], lscombos_sel_pilot$pleiorate[j], lscombos_sel_pilot$delmu[j], as.integer(round(lscombos_sel_pilot$nloci[j])), lscombos_sel_pilot$locisigma[j], as.integer(round(lscombos_sel_pilot$delchr[j])), lscombos_sel_pilot$wsd[j], j, intern=T))
  }
stopCluster(cl)
