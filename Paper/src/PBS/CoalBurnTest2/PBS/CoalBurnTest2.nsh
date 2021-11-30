#!/sw7/RCC/NimrodG/embedded-1.9.0/bin/nimexec

# Modified from original script by David.Green@uq.edu.au
# More information at: https://github.com/UQ-RCC/nimrod-embedded
#
# Submit this script as a PBS job script using the PBSPro qsub command.
# The nimexec command will parse it into what is required.

#Fix up Account String
#PBS -A qris-uq
#
#Use more resources if you need.
#PBS -l select=4:ncpus=24:mem=120GB:ompthreads=1 
#40 nodes, 24 cores, 120GB per node
#PBS -l walltime=48:00:00
#
#Job name for ease of recognition
#PBS -N CoalBurnTest2
#

# There are additional directives for Nimrod to interpret with #NIM at the start of each line.
# Tell Nimrod to use this as the shell for the job proper when it has parsed this file.
#NIM shebang /bin/bash

# =============================================================================
# Tell Nimrod what range of parameter values you want to use "1 to 1638 step 1" 
# =============================================================================

#The parameters for the 100 latin squares are rows in the input file.
#NIM parameter LS integer range from 1 to 64 step 1

#Repeat 50 times for each hypercube with a different SEED value
#NIM parameter SEED integer range from 1 to 10 step 1


# Just checking that something did not go wrong with assignment of the J values.
if [ -z "${NIMROD_VAR_LS}" ]; then
        echo "\$NIMROD_VAR_LS isn't set, cannot continue..." 
        exit 2
fi

if [ -z "${NIMROD_VAR_SEED}" ]; then
        echo "\$NIMROD_VAR_SEED isn't set, cannot continue..." 
        exit 2
fi

#Where you submit this job from will be the value of $PBS_O_WORKDIR
echo "PBS_O_WORKDIR is ${PBS_O_WORKDIR}" 
#Everything you need should be located relative to PBS_O_WORKDIR, or else a full path
#Set the cd to TMPDIR for writing SLiM output
cd ${TMPDIR}

#=====================
#Modify these to suit.
#=====================

# Always run the entire parameter range cause nimrod can do them in any order.
# See the -f test below about skipping the ones we have already done.
RUNNAME="CoalBurnTest2_singleshot" 

OUTFILE="${PBS_O_WORKDIR}/Outputs/TEST_${NIMROD_VAR_SEED}_${NIMROD_VAR_LS}.txt" 
echo "${OUTFILE}" 

if [ -f ${OUTFILE} ]; then
  echo "Output file ${OUTFILE} already exists. Skipping this index value ${NIMROD_VAR_SEED} ${NIMROD_VAR_LS}" 
  exit 0
fi


RSCRIPTNAME="${PBS_O_WORKDIR}/R/${RUNNAME}.R" 

module purge
module load R/4.0.0

Rscript $RSCRIPTNAME ${NIMROD_VAR_LS} ${NIMROD_VAR_SEED} 
