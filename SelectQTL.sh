#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p ghpc_v3                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 3                       # Number of CPU cores
#SBATCH --mem=80024                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J SelQTL_job            # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.err    # STDERR
#SBATCH -t 48:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #

NumQTL=${1?Error: QTL number not given}
NumNeu=${2?Error: Neutral marker number not given}
NumRep=${3?Error: replicate number not given}

wd=$(pwd)

for i in $(seq 1 $NumRep)
do
Rscript SelectQTLBeforeSimulation.R $NumQTL $NumNeu $i

done
