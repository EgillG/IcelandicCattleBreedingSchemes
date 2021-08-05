#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #    
#--------------------------------------------------------------------------#
#SBATCH -p ghpc_v3                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 2                       # Number of CPU cores
#SBATCH --mem=20000                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J ocs_sim_job            # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.err    # STDERR
#SBATCH -t 144:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #

replicate=${1?Error: replicate number not given}
method=${2?Error: method not given}
NumBulls=${3?Error: Number of bulls not given}

mkdir evaSim
mkdir $method

Rscript IceSim_Test_04_08_2021.R $replicate $method $NumBulls

