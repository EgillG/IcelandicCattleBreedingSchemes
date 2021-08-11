#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #    
#--------------------------------------------------------------------------#
#SBATCH -p ghpc_v1                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 1                       # Number of CPU cores
#SBATCH --mem=1024                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J template_job            # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.err    # STDERR
#SBATCH -t 1:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM

# wd is working drive
wd=$(pwd)
replicate=${1?Error: replicate number not given}
method=${2?Error: method not given}
NumBulls=${3?Error: number of bulls not given}
echo $replicate $method $NumBulls

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

if [ $method == "M1" ] || [ $method == "M2" ] || [ $method == "M1D" ] || [ $method == "M2D" ] || [ $method == "M1D_Est" ] || [ $method == "M2D_Est" ] || [ $method == "M1R" ] || [ $method == "M2R" ]
then
echo method $method starting GMATRIX
python3 gmatPar.py $wd $method

echo " Gmatrix Job started at $(date '+%y-%m-%d %H:%M:%S')"

JOBNAME=IceSim
PARFILE=$SLURM_SUBMIT_DIR/Gmatrix/gmat.par
PROG=/usr/home/qgg/gs/public/invg-md-all-v8/invgmatrix

cp $PROG $TMPDIR/invgmatrix
cp $PARFILE $TMPDIR/par.dat

cd $TMPDIR/
echo running GMATRIX
ulimit -s unlimited
echo $SLURM_SUBMIT_DIR
./invgmatrix

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID

echo "Gmatrix Job completed at $(date '+%y-%m-%d %H:%M:%S')"

elif [ $method == "H1" ] || [ $method == "H2" ] || [ method == "H3" ]
then
# here method is an input to the haplotypes program, to tell it how long haplotypes to use
python3 Haplotpes_MoBOPS_8_15_02_2021.py $method

echo " Gmatrix Job started at $(date '+%y-%m-%d %H:%M:%S')"

JOBNAME=IceSim
PARFILE=$SLURM_SUBMIT_DIR/Gmatrix/gmat.par
PROG=/usr/home/qgg/gs/public/invg-md-all-v8/invgmatrix

cp $PROG $TMPDIR/invgmatrix
cp $PARFILE $TMPDIR/par.dat

cd $TMPDIR/
echo running GMATRIX
ulimit -s unlimited
echo $SLURM_SUBMIT_DIR
./invgmatrix

cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID

echo "Gmatrix Job completed at $(date '+%y-%m-%d %H:%M:%S')"

fi

# This program makes the EVA.prm file
python3 EVApar.py $method $NumBulls $wd

prg=/opt/ghpc/eva/bin/eva
$prg EVA.prm