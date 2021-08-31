#!/bin/bash -l

cd $PBS_JOBFS
SECONDS=0

R --file=/home/$USER/SLiM/Scripts/Tests/GeneticConstraints/Sept2021_GadiRun/R/constraints_singleshot_sublaunch.R

cat ./
cat /$TMPDIR/out_stabsel_pos.csv >> /30days/$USER/out_stabsel_pos.csv


DURATION=$SECONDS

echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
