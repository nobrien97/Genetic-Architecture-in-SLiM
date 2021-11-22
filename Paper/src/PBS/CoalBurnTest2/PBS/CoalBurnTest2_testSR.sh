#!/bin/bash -l

module load R/4.0.0
module load python3/3.9.2

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2

RSCRIPTNAME=$HOME/tests/CoalBurnTest2/R/CoalBurnTest2_singleshot.R

echo "Running modelindex = $MODELINDEX, seedindex = $SEED...\n"
Rscript ${RSCRIPTNAME} ${MODELINDEX} ${SEED}


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seedindex = $SEED, finished successfully!\n"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."


