 ##############################################################################################################
#  Coalescent (backwards time) burn-in test
##############################################################################################################

#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.
#  Thanks to David Green at the UQ RCC


# Parse parameters: seed and latin square row number


args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 2 ) {
  cat("Need 2 command line parameters: LS_COMBO SEED\n")
  q()
}

# combo and seed
i <- as.numeric(args[1])
j <- as.numeric(args[2])

# Env variables
HOME <- Sys.getenv("HOME")

# Seeds generated with seedgenerator 
seeds <- read.csv(paste0(HOME, "/tests/CoalBurnTests/R/seeds_coalburntest.csv"), header = T)
combos <- read.csv(paste0(HOME, "/tests/CoalBurnTests/R/lhc_coalburntest.csv"), header = T)


# Check if we are a control or tree model
if ( combos$tree[j] ) {
  # Sublaunch Python to generate the tree sequence (conda sucks)
  system(sprintf("$HOME/.conda/envs/msslim/bin/python3.9 $HOME/tests/CoalBurnTests/R/burnin.py %i %s", 
                j, 
                as.character(seeds$Seed[i]), 
                intern = T))

  # Run SLiM
  slim_out <- system(sprintf("$HOME/slim/slim -s %s -d Ne=%i -d rwide=%f -d modelindex=%i $HOME/tests/CoalBurnTests/slim/polygen_maint_coalload.slim", 
                              as.character(seeds$Seed[i]), 
                              as.integer(ceiling(combos[j,]$Ne)), 
                              combos[j,]$rwide, 
                              j, intern=T))
} else {
# Sublaunch SLiM with the appropriate values
slim_out <- system(sprintf("$HOME/slim/slim -s %s -d Ne=%i -d rwide=%f -d modelindex=%i $HOME/tests/CoalBurnTests/slim/polygen_maint.slim", 
                           as.character(seeds$Seed[i]), 
                           as.integer(round(combos[j,]$Ne)), 
                           combos[j,]$rwide, 
                           j, intern=T))

}


