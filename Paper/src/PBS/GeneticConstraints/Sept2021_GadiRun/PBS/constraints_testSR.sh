#!/bin/bash -l

module load R/4.0.0

cd $PBS_JOBFS
SECONDS=0

# Rename the first and second arguments passed to this single shot script for clarity 
MODELINDEX=$1
SEED=$2

RSCRIPTNAME=$HOME/SLiM/Tests/Sept2021_GadiRun/R/constraints_singleshot_sublaunch.R

echo "Running modelindex = $MODELINDEX, seedindex = $SEED...\n"
Rscript ${RSCRIPTNAME} ${MODELINDEX} ${SEED}

echo "Copying output..."
cat ./out_stabsel_pos* >> /scratch/ht96/nb9894/slim_SeptGadiTest/out_stabsel_pos.csv
cat ./out_stabsel_burnin* >> /scratch/ht96/nb9894/slim_SeptGadiTest/out_stabsel_burnin.csv
cat ./out_stabsel_means* >> /scratch/ht96/nb9894/slim_SeptGadiTest/out_stabsel_means.csv
cat ./out_stabsel_opt* >> /scratch/ht96/nb9894/slim_SeptGadiTest/out_stabsel_opt.csv
cat ./out_stabsel_muts* >> /scratch/ht96/nb9894/slim_SeptGadiTest/out_stabsel_muts.csv
cat ./out_stabsel_dict* >> /scratch/ht96/nb9894/slim_SeptGadiTest/out_stabsel_dict.csv


DURATION=$SECONDS
echo "Run modelindex = $MODELINDEX, seedindex = $SEED, finished successfully!\n"
echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."


