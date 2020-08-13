##############################################################################################################
#  Pilot run stabilising selection model script with arguments coming from latin square, for use on cluster
##############################################################################################################

#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

#NEED TO PROCESS 2 PARAMETERS PASSED IN  REPEAT(i.e. SEED) and LATIN SQUARE ROW NUMBER
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 2 ) {
  cat("Need 2 command line parameters i.e. SEED LS_COMBO\n")
  q()
}

row_seed     <- as.numeric(args[1])
row_combo    <- as.numeric(args[2])

# Environment variables

USER <- Sys.getenv('USER')


# Load LHS samples - For this test run, we are doing 8 LHC samples and 3 seeds

lscombos_null_pilot <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Pilot/AIM1/Inputs/lscombos_null_pilot.csv"), header = T)

rsample <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Pilot/AIM1/Inputs/seeds_pilot.csv"), header = T) # .csv in a single column

#Run SLiM, defining parameter sets according to LHC samples in command line

# 8 Parameters for Aim 1 and 2, 9 for Aim 3 (wsd):
# Ne, rwide, pleio_cov, pleiorate, delmu, nloci, locisigma, delchr, wsd

i <- row_seed
j <- row_combo

slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d nloci=%i -d locisigma=%f -d delchr=%i -d modelindex=%i /home/$USER/SLiM/Scripts/Nimrod/Pilot/AIM1/Inputs/null8T100L.slim",
                       as.character(rsample$Seed[i]), as.integer(round(lscombos_null_pilot$Ne[j])), as.character(lscombos_null_pilot$rwide[j]), lscombos_null_pilot$pleiocov[j], lscombos_null_pilot$pleiorate[j], lscombos_null_pilot$delmu[j], as.integer(round(lscombos_null_pilot$nloci[j])), lscombos_null_pilot$locisigma[j], as.integer(round(lscombos_null_pilot$delchr[j])), j, intern=T))
