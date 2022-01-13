#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p ghpc_v3                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 1                       # Number of CPU cores
#SBATCH --mem=80024                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J template_job            # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.err    # STDERR
#SBATCH -t 48:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #

NumRep=${1?Error: replicate number not given}
wd=$(pwd)

for i in $(seq 1 $NumRep)
do
python3 QMSimtoPlink_clusterVersion.py $wd $i
if [ $i -lt 10 ]
then
awk 'NR>1{print $2,substr($1,2,5),$3,0}' r_IceSim/lm_mrk_00${i}.txt > data_${i}.map
plink-1.9-rc --file data_${i} --cow --freq --out plink_${i}
awk 'NR>1{print $2}' plink_${i}.frq | shuf -n 60000 > extr
plink-1.9-rc --file data_${i} --cow --extract extr --recode --out data_${i}_extr
awk '{print $2,NR}' data_${i}_extr.map > remap

join -1 1 -2 2 <(sort extr) <(sort -k2,2 plink_${i}.frq) | tr ' ' '\t' | sort -k2,2n -k1,1n |awk '{print NR,$2,$3,$4,$5,$6}' | tr ' ' '\t' > temp
mv temp plink_${i}.frq
mv data_${i}_extr.map data_${i}.map
mv data_${i}_extr.ped data_${i}.ped
cat data_${i}.map | awk '{print $1,NR,$3,$4}' | tr ' ' '\t' > temp
mv temp data_${i}.map

else
awk 'NR>1{print $2,substr($1,2,5),$3,0}' r_IceSim/lm_mrk_0${i}.txt > data_${i}.map
plink-1.9-rc --file data_${i} --cow --freq --out plink_${i}
awk 'NR>1{print $2}' plink_${i}.frq | shuf -n 60000 > extr
plink-1.9-rc --file data_${i} --cow --extract extr --recode --out data_${i}_extr
awk '{print $2,NR}' data_${i}_extr.map > remap
plink-1.9-rc --file data_${i}_extr --cow --update-map remap --recode --out temp
mv temp.map data_${i}_extr.map
mv temp.ped data_${i}_extr.ped
rm temp.map temp.ped

join -1 1 -2 2 <(sort extr) <(sort -k2,2 plink_${i}.frq) | tr ' ' '\t' | sort -k2,2n -k1,1n |awk '{print NR,$2,$3,$4,$5,$6}' | tr ' ' '\t' > temp
mv temp plink_${i}.frq
mv data_${i}_extr.map data_${i}.map
mv data_${i}_extr.ped data_${i}.ped
cat data_${i}.map | awk '{print $1,NR,$3,$4}' | tr ' ' '\t' > temp
mv temp data_${i}.map

fi
done
