

# Tau value testing + time trial for selection model


TMPDIR <- Sys.getenv('TMPDIR')
USER <- Sys.getenv('USER')
PBSIND <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

library(foreach)
library(doParallel)
library(future)

rsample <- c(78478, 167, 244978)

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

if (PBSIND == 1) {
	foreach(i=rsample) %dopar%
		system.time(system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=8000 -d locisigma=10.0 -d pleiorate=0.5 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=0.5 -d tau=100.0 /home/$USER/SLiM/Scripts/tautest/stabsel_recom_8T100L_fitness.slim",
									as.character(i), intern=T)))
} else if (PBSIND == 2) {
	foreach(i=rsample) %dopar%
		system.time(system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=8000 -d locisigma=10.0 -d pleiorate=0.5 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=0.5 -d tau=1000.0 /home/$USER/SLiM/Scripts/tautest/stabsel_recom_8T100L_fitness.slim",
									as.character(i), intern=T)))
} else {
	foreach(i=rsample) %dopar%
		system.time(system(sprintf("/home/$USER/SLiM/slim -s %s -d Ne=8000 -d locisigma=10.0 -d pleiorate=0.5 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=0.5 -d tau=10000.0 /home/$USER/SLiM/Scripts/tautest/stabsel_recom_8T100L_fitness.slim",
									as.character(i), intern=T)))  
    }
	
stopCluster(cl)
