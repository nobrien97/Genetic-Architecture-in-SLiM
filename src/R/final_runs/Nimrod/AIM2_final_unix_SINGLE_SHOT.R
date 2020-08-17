##########################################################################################################################
#  Final run drift with independent assortment model script with arguments coming from latin square, for use on cluster  #
##########################################################################################################################

#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Thanks to David Green (David.Green@uq.edu.au) for getting this to work on Nimrod

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


# Load LHS samples - Production run is 100 samples, 50 seeds

lscombos_recom <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Final/AIM2/Inputs/lscombos_null.csv"), header = T)

rsample <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Nimrod/Final/AIM2/Inputs/seeds.csv"), header = T) # .csv in a single column

#Run SLiM, defining parameter sets according to LHC samples in command line

# 8 Parameters for Aim 1 and 2, 9 for Aim 3 (wsd):
# Ne, rwide, pleio_cov, pleiorate, delmu, nloci, locisigma, delchr, wsd

i <- row_seed
j <- row_combo

slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%s -d pleio_cov=%f -d pleiorate=%f -d delmu=%f -d nloci=%i -d locisigma=%f -d delchr=%i -d modelindex=%i /home/$USER/SLiM/Scripts/Nimrod/Final/AIM2/Inputs/null_recom_8T100L.slim",
                       as.character(rsample$Seed[i]), as.integer(round(lscombos_recom$Ne[j])), as.character(lscombos_recom$rwide[j]), lscombos_recom$pleiocov[j], lscombos_recom$pleiorate[j], lscombos_recom$delmu[j], as.integer(round(lscombos_recom$nloci[j])), lscombos_recom$locisigma[j], as.integer(round(lscombos_recom$delchr[j])), j, intern=T))
