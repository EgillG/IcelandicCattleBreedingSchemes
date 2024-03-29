#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #    
#--------------------------------------------------------------------------#
#SBATCH -p ghpc_v2                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 4                       # Number of CPU cores
#SBATCH --mem=80000                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J ocs_sim_job            # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.err    # STDERR
#SBATCH -t 240:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #

replicate=${1?Error: replicate number not given}
method=${2?Error: method not given}
NumBulls=${3?Error: Number of bulls not given}

mkdir ${method}_${replicate}

cp IceSim_main.R EVApar.py EVA.sh gmatPar.py haplo.py G_matrixPreparation_3.py r_dmu5 r_dmu4 DMU.sh dmuSS.DIR dmuPBLUP.DIR dmuGBLUP.DIR ${method}_${replicate}
cd ${method}_${replicate}
mkdir Gmatrix
Rscript IceSim_main.R $replicate $method $NumBulls

