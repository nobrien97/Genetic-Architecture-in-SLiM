#!/bin/bash -l
#PBS -q normal
#PBS -N slim_gadiTest
#PBS -l walltime=05:00:00
#PBS -l ncpus=624
#PBS -l mem=768GB
#PBS -l jobfs=1000GB
#PBS -l storage=scratch/ht96+gdata/ht96

# Make sure we're at the right place so we can find the bash script to run
cd $PBS_O_WORKDIR

# Make output folder
mkdir /scratch/ht96/nb9894/CoalBurnTest2

# Analogous to UQ Tinaroo embedded Nimrod
# Use 1 core per SLiM run
module load nci-parallel/1.0.0
export ncores_per_task=1
export ncores_per_numanode=12
CMDS_PATH=$HOME/tests/CoalBurnTest2/PBS/cmds.txt

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file $CMDS_PATH --timeout 86400

# Combine output into a single file
cd /scratch/ht96/nb9894/CoalBurnTest2/
cat ./out_stabsel_pos* >> ./out_stabsel_pos.csv
cat ./out_stabsel_burnin* >> ./out_stabsel_burnin.csv
cat ./out_stabsel_means* >> ./out_stabsel_means.csv
cat ./out_stabsel_opt* >> ./out_stabsel_opt.csv
cat ./out_stabsel_muts* >> ./out_stabsel_muts.csv
cat ./out_stabsel_dict* >> ./out_stabsel_dict.csv
cat ./out_stabsel_time* >> ./out_stabsel_time.csv
cat ./tree_time_* >> ./tree_time.csv

# Zip LD matrices
/bin/zip -q -Z bzip2 ./out_stabsel_ld.zip ./out_stabsel_ld* 

# Delete loose files with seed and model indices
find -regex ".*[0-9]*_*[0-9].csv+" -delete

