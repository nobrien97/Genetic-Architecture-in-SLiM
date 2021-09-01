#!/bin/bash -l
#PBS -q normal
#PBS -N slim_gadiTest
#PBS -l walltime=10:00:00
#PBS -l ncpus=2880
#PBS -l mem=1140GB
#PBS -l jobfs=1000GB
#PBS -l storage=scratch/ht96+gdata/ht96

# Make sure we're at the right place so we can find the bash script to run
cd $PBS_O_WORKDIR

# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0
export ncores_per_task=1
export ncores_per_numanode=12
CMDS_PATH=$HOME/SLiM/Tests/Sept2021_GadiRun/PBS/cmds.txt

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file $CMDS_PATH --timeout 86400

# Copy the output files to project data, and combine them into a single file
cd $PBS_JOBFS
mkdir /scratch/ht96/slim_SeptGadiTest
cat ./out_stabsel_pos* >> /scratch/ht96/slim_SeptGadiTest/out_stabsel_pos.csv
cat ./out_stabsel_burnin* >> /scratch/ht96/slim_SeptGadiTest/out_stabsel_burnin.csv
cat ./out_stabsel_means* >> /scratch/ht96/slim_SeptGadiTest/out_stabsel_means.csv
cat ./out_stabsel_opt* >> /scratch/ht96/slim_SeptGadiTest/out_stabsel_opt.csv
cat ./out_stabsel_muts* >> /scratch/ht96/slim_SeptGadiTest/out_stabsel_muts.csv
cat ./out_stabsel_dict* >> /scratch/ht96/slim_SeptGadiTest/out_stabsel_dict.csv

