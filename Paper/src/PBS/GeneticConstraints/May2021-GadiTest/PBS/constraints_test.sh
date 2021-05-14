#!/bin/bash -l
#PBS -q normal
#PBS -p $PROJECT
#PBS -N slim_gadiTest
#PBS -l walltime=220:00:00
#PBS -l ncpus=15360
#PBS -l mem=61440GB
#PBS -l jobfs=128000GB
#PBS -l storage=scratch/ht96+gdata/ht96

cd $TMPDIR

module load nci-parallel/1.0.0
export ncores_per_task=1
export ncores_per_numanode=12

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file cmds.txt --timeout 86400