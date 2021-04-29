##############################################################################################################
#  Genetic Constraint run: 48 replicates, stabilising selection, 40 genes, 30 QTLs, recom, locisigma, gen con
##############################################################################################################


#  Parallel script modified from SLiM-Extras example R script, info at
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.

# Environment variables

USER <- Sys.getenv('USER')
ARR_INDEX <- as.numeric(Sys.getenv('PBS_ARRAY_INDEX'))


# Parallelisation libraries 

library(foreach)
library(doParallel)
library(future)


seeds <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Tests/GeneticConstraints/Apr2021-LS_LMHgc_r/R/seeds_48.csv"), header = T)
combos <- read.csv(paste0("/home/",USER,"/SLiM/Scripts/Tests/GeneticConstraints/Apr2021-LS_LMHgc_r/R/combos.csv"), header = T)

# Set which runs to do according to node

switch (object,
  case = action
)

switch (ARR_INDEX,
  { combos <- combos[1:5,] },
  { combos <- combos[6:9,] },
  { combos <- combos[10:13,] },
  { combos <- combos[14:18,] }
)


cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

#Run SLiM
foreach(i=1:nrow(combos)) %:%
  foreach(j=seeds$Seed) %dopar% {
	# Use string manipulation functions to configure the command line args, feeding from a data frame of seeds
	# then run SLiM with system(),
    	slim_out <- system(sprintf("/home/$USER/SLiM/slim -s %s -d locisigma=%f -d rwide=%f -d con_props=%s -d modelindex=%i /home/$USER/SLiM/Scripts/Tests/GeneticConstraints/Apr2021-LS_LMHgc_r/slim/polygen_maint.slim", as.character(j), combos[i,]$locisigma, combos[i,]$rwide, combos[i,]$gencon, combos[i,]$modelindex, intern=T))
  }
stopCluster(cl)

