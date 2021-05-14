 ##############################################################################################################
#  Gadi test run, 6 LHC parameters
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Parse parameters: seed and latin square row number
 
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 2 ) {
  cat("Need 2 command line parameters: SEED LS_COMBO\n")
  q()
}

row_seed <- as.numeric(args[1])
row_combo <- as.numeric(args[2])
 
# Environment variables
USER <- Sys.getenv('USER')


# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)

# Seeds generated with seedgenerator, from 
seeds <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Tests/May2021-GadiTest/R/seeds_gaditest.csv"), header = T)
combos <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Tests/May2021-GadiTest/R/lhc_samples"), header = T)

i <- row_seed
j <- row_combo

slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=%i -d rwide=%f -d nloci=%i -d locisigma=%f -d con_props=%s -d width=%f -d opt=%f -d modelindex=%i /home/$USER/SLiM/Scripts/Tests/May2021-GadiTest/slim/polygen_maint.slim", 
                           as.character(i), 
                           as.integer(round(combos[j,]$Ne)), 
                           combos[j,]$rwide, 
                           as.integer(round(combos[j,]$nloci)), 
                           combos[j,]$locisigma, 
                           combos[j,]$K, 
                           combos[j,]$width, 
                           combos[j,]$opt, 
                           j, intern=T))


