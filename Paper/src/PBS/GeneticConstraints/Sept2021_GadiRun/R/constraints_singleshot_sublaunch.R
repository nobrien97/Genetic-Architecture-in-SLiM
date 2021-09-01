 ##############################################################################################################
#  Gadi test run, 6 LHC parameters
##############################################################################################################

#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Parse parameters: seed and latin square row number
 
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 2 ) {
  cat("Need 2 command line parameters: LS_COMBO SEED\n")
  q()
}

row_combo <- as.numeric(args[1])
row_seed <- as.numeric(args[2])

# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)

# Seeds generated with seedgenerator 
seeds <- read.csv("$HOME/SLiM/Tests/Sept2021_GadiRun/R/seeds_gaditest.csv", header = T)
combos <- read.csv("$HOME/SLiM/Tests/Sept2021_GadiRun/R/lhc_gaditest.csv", header = T)

i <- row_seed
j <- row_combo

# Sublaunch SLiM with the appropriate values
slim_out <- system(sprintf("$HOME/SLiM/slim -s %s -d Ne=%i -d rwide=%f -d nloci=%i -d locisigma=%f -d con_props=%s -d width=%f -d optShift=%f -d modelindex=%i $HOME/SLiM/Tests/Sept2021_GadiRun/slim/polygen_maint.slim", 
                           as.character(i), 
                           as.integer(round(combos[j,]$Ne)), 
                           combos[j,]$rwide, 
                           as.integer(round(combos[j,]$nloci)), 
                           combos[j,]$locisigma, 
                           combos[j,]$K, 
                           combos[j,]$width, 
                           combos[j,]$optShift, 
                           j, intern=T))


