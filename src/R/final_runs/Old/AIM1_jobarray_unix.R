##############################################################################################################
#  Final Null model script with arguments coming from latin square, for use on cluster
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


# Load LHS samples - We are splitting into 100 parts, so each node does 1 row of the lscombos_null file (100 total runs/node)

lscombos_null <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/lscombos_null.csv"), header = T)

rsample <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/seeds.csv"), header = F)
rsample <- as.character(as.vector(t(rsample)))

# - Multinode parallelism: Use the array index to tell us what range in lscombos_null to evaluate

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM, defining parameter sets according to LHS in command line - j is from (PBSIND*5 )-4 to PBSIND*5


# 8 Parameters for Aim 1 and 2, 9 for Aim 3 (wsd):
# Ne, rwide, pleio_cov, pleiorate, delmu, nloci, locisigma, delchr


foreach(i=rsample) %:%
  foreach(j=PBSIND) %dopar% {
    # Use string manipulation functions to configure the command line args, feeding from a data frame of latin square combinations
    # then run SLiM with system(),
    slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d nloci=%i -d locisigma=%f -d delchr=%i -d modelindex=%i /home/$USER/SLiM/Scripts/null8T100L.slim", 
                               as.character(i), as.integer(round(lscombos_null$Ne[j])), as.character(lscombos_null$rwide[j]), lscombos_null$pleiocov[j], lscombos_null$pleiorate[j], lscombos_null$delmu[j], as.integer(round(lscombos_null$nloci[j])), lscombos_null$locisigma[j], as.integer(round(lscombos_null$delchr[j])), j, intern=T))
  }
stopCluster(cl)
