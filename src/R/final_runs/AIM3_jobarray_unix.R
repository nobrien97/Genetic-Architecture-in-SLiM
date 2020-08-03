##############################################################################################################
#  Final stabilising selection model script with arguments coming from latin square, for use on cluster
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Get environment variables

TMPDIR <- Sys.getenv('TMPDIR')
USER <- Sys.getenv('USER')
PBSIND <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

print(PBSIND)


# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


# Load LHS samples - We are splitting into 100 parts, so each node does 1 row of the lscombos_sel file (100 total runs/node)

lscombos_sel <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/lscombos_sel.csv"), header = T)

rsample <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/seeds.csv"), header = F)
rsample <- as.character(as.vector(t(rsample)))

# - Multinode parallelism: Use the array index to tell us what range in lscombos_sel to evaluate

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to LHS in command line - j is the PBS_ARRAY_INDEX value. Each node will do one row (100 seeds)


# 8 Parameters for Aim 1 and 2, 9 for Aim 3 (wsd):
# Ne, rwide, pleio_cov, pleiorate, delmu, nloci, locisigma, delchr, wsd


foreach(i=rsample) %:%
  foreach(j=PBSIND) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d nloci=%i -d locisigma=%f -d delchr=%i -d wsd=%f -d modelindex=%i /home/$USER/SLiM/Scripts/stabsel_recom_8T100L.slim", 
                               as.character(i), as.integer(round(lscombos_sel$Ne[j])), as.character(lscombos_sel$rwide[j]), lscombos_sel$pleiocov[j], lscombos_sel$pleiorate[j], lscombos_sel$delmu[j], as.integer(round(lscombos_sel$nloci[j])), lscombos_sel$locisigma[j], as.integer(round(lscombos_sel$delchr[j])), lscombos_sel$wsd[j], j, intern=T))
  }
stopCluster(cl)
