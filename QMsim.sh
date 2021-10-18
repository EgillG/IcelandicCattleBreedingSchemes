#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p ghpc_v3                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 2                       # Number of CPU cores
#SBATCH --mem=10024                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J template_job            # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.err    # STDERR
#SBATCH -t 4:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #

NumRep=${1?Error: replicate number not given}
NumChrom=${2?Error: Number of chromosomes not given}
NumMark=${3?Error: Number of markers per chromosome not given}
SizeChrom=${4?Error: Size of chromosomes not given. Supply size in centiMorgans}
wd=$(pwd)
rm r_IceSim/*

python3 QMsimPar.py $NumRep $NumChrom $NumMark $SizeChrom $wd

/opt/ghpc/QMSim_Linux/QMSim IceSim.prm -o

wait
#while [[ ! -f /r_IceSim/stat_rep.txt ]]
#do
#sleep 1
#done

for i in $(seq 1 $NumRep)
do
python3 QMSimtoPlink_clusterVersion.py $wd $i
if [ $i -lt 10 ]
then
awk 'NR>1{print $2,substr($1,2,5),$3,0}' r_IceSim/lm_mrk_00${i}.txt > data_${i}.map
plink-1.9-rc --file data_${i} --cow --freq --out plink_${i}
else
awk 'NR>1{print $2,substr($1,2,5),$3,0}' r_IceSim/lm_mrk_0${i}.txt > data_${i}.map
plink-1.9-rc --file data_${i} --cow --freq --out plink_${i}
fi
done


rm plink_*.nosex plink_*.log
