rep=${1?Error: Replicate number not given}
num=${2?Error: Number of bulls not given}

#for i in $(seq 1 $rep)
#do
i=$rep
#sbatch ocs.sh $i Ped $num
#sbatch ocs.sh $i M1 $num
sbatch ocs.sh $i M2 $num
#sbatch ocs.sh $i M1D $num
sbatch ocs.sh $i M2D $num
#sbatch ocs.sh $i M1R $num
#sbatch ocs.sh $i M2R $num
#sbatch ocs.sh $i M1_05 $num
#done
