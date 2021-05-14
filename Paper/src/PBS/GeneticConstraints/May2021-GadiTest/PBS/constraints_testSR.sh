#!/bin/bash -l

cd $TMPDIR
SECONDS=0

R --file=/home/$USER/SLiM/Scripts/Tests/GeneticConstraints/Apr2021-LS_LMHgc_r/R/constraints_sublaunch.R

cat /$TMPDIR/out_stabsel_means.csv >> /30days/$USER/out_stabsel_means.csv

cat /$TMPDIR/out_stabsel_muts.csv >> /30days/$USER/out_stabsel_muts.csv

cat /$TMPDIR/out_stabsel_burnin.csv >> /30days/$USER/out_stabsel_burnin.csv

cat /$TMPDIR/out_stabsel_opt.csv >> /30days/$USER/out_stabsel_opt.csv

cat /$TMPDIR/out_stabsel_dict.csv >> /30days/$USER/out_stabsel_dict.csv

cat /$TMPDIR/out_stabsel_pos.csv >> /30days/$USER/out_stabsel_pos.csv


DURATION=$SECONDS

echo "$(($DURATION / 3600)) hours, $((($DURATION / 60) % 60)) minutes, and $(($DURATION % 60)) seconds elapsed."
