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

# combo and seed argument from command line
i <- as.numeric(args[1])
j <- as.numeric(args[2])

# Env variables
HOME <- Sys.getenv("HOME")

# Seeds generated with seedgenerator 
seeds <- read.csv(paste0(HOME, "/tests/CoalBurnTests/R/seeds_coalburntest.csv"), header = T)
combos <- read.csv(paste0(HOME, "/tests/CoalBurnTests/R/lhc_coalburntest.csv"), header = T)

# Set timer
time <- system.time(
  # Check if we are a control or tree model
  if(combos$tree[i]){
    # Sublaunch Python to generate the tree sequence (conda sucks)
    system(sprintf("$HOME/.conda/envs/msslim/bin/python3.9 $HOME/tests/CoalBurnTests/R/burnin.py %i %s",
                  i,
                  as.character(seeds$Seed[j]),
                  intern=T))

    # Run SLiM
    slim_out <- system(sprintf("$HOME/slim/slim -s %s -d Ne=%i -d rwide=%f -d modelindex=%i $HOME/tests/CoalBurnTests/slim/polygen_maint_coalload.slim",
                                as.character(seeds$Seed[j]),
                                as.integer(ceiling(combos[i,]$Ne)),
                                combos[i,]$rwide,
                                i,intern=T))
  } else{
    # Sublaunch SLiM with the appropriate values
    slim_out <- system(sprintf("$HOME/slim/slim -s %s -d Ne=%i -d rwide=%f -d modelindex=%i $HOME/tests/CoalBurnTests/slim/polygen_maint.slim",
                            as.character(seeds$Seed[j]),
                            as.integer(round(combos[i,]$Ne)),
                            combos[i,]$rwide,
                            i,intern=T))

  }
)[3] # Show elapsed time only

# Save the total time taken to run the simulation
write.table(time, paste0("/scratch/user/uqnobri4/CoalBurnTests/time_", seeds$Seed[j], "_", i), 
          row.names = F, col.names = F)

