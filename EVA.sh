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
echo submit directory $SLURM_SUBMIT_DIR

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR

if [ $method == "M1" ] || [ $method == "M2" ] || [ $method == "M1D" ] || [ $method == "M2D" ] || [ $method == "M1D_Est" ] || [ $method == "M2D_Est" ] || [ $method == "M1R" ] || [ $method == "M2R" ]
then
echo method $method starting GMATRIX
python3 gmatPar.py $wd $method Yes
#awk '{ if (substr($2,2,1)!= "e") print $0}' mapbase.dat > temp
#mv temp mapbase.dat

mv gmat.dat gmat.id gmat.map mapbase.dat $SLURM_SUBMIT_DIR/${method}_${replicate}/Gmatrix/

echo " Gmatrix Job started at $(date '+%y-%m-%d %H:%M:%S')"

JOBNAME=IceSim
PARFILE=$SLURM_SUBMIT_DIR/${method}_${replicate}/Gmatrix/gmat.par
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

if [ $method == "H1" ] || [ $method == "H2" ] || [ $method == "H3" ] || [ $method == "H4" ]
then
echo Haplotype based matrices
tail -n +2 Gmatrix/gmat > temp
mv temp Gmatrix/gmat
fi
cd $SLURM_SUBMIT_DIR

mv ${method}_${replicate}/Gmatrix/gmat ${method}_${replicate}/Gmatrix.gmat
# This program makes the EVA.prm file
python3 EVApar.py $method $NumBulls $wd
cd ${method}_${replicate}
#cut -f1,4,6 evaIn.txt > SelCands
# Use grep to remove animals that are not candidates for selection.
if [ $method != "Ped" ]
then
grep -wf  <(cut -f1 SelCands) Gmatrix.gmat > ReducedMatrix
fi
echo $wd
python3 G_matrixPreparation_3.py
mkdir evaSim

echo "EVA Job started at $(date '+%y-%m-%d %H:%M:%S')"

prg=/opt/ghpc/eva/bin/eva
$prg EVA.prm
echo lines in Gmatrix are:
wc -l Gmatrix.gmat
rm Gmatrix.gmat ReducedMatrix
rm haplomatrix*

echo "EVA Job completed at $(date '+%y-%m-%d %H:%M:%S')"
rm GINV*
rm Gmatrix/gmat.dat
rm Gmatrix/invgmat
rm Gmatrix/gmat
