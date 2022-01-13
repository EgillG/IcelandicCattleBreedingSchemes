method=${1?Error: method not given to DMU.sh}
singleS=${2?Error: Single step not given to DMU.sh}
replicate=${3?Error: replicate not given to DMU.sh}
echo " Gmatrix Job started at $(date '+%y-%m-%d %H:%M:%S')"
echo running GMATRIX to prepare for single step breeding value estimation with DMU.

TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR


#PARFILE=$SLURM_SUBMIT_DIR/${method}/Gmatrix/gmat.par
PARFILE=Gmatrix/gmat.par
PROG=/usr/home/qgg/gs/public/invg-md-all-v8/invgmatrix

cp $PROG $TMPDIR/invgmatrix
cp $PARFILE $TMPDIR/par.dat

cd $TMPDIR/
echo running GMATRIX
ulimit -s unlimited
echo $SLURM_SUBMIT_DIR
./invgmatrix

cd $SLURM_SUBMIT_DIR/${method}_${replicate}/
rm -rf /scratch/$USER/$SLURM_JOBID

echo "Gmatrix Job completed at $(date '+%y-%m-%d %H:%M:%S')"
cut -f 1 -d ' ' Gmatrix/gmat.dat > genotyped
if [ $singleS = 'Yes' ]
then
echo running single step GBLUP in DMU
bash r_dmu5 dmuSS
awk '{print $5,$8}' dmuSS.SOL > temp
mv temp dmuSS.SOL
elif  [ $singleS = 'No' ]
then
echo running GBLUP in DMU
#join -1 1 -2 1 -o 1.1,2.2,2.3 <(cut -f 1 -d ' ' Gmatrix/gmat.dat | sort -k1,1 ) <(cat phenotypes.txt | sort -k1,1) > tempPhen
#mv tempPhen phenotypes.txt
bash r_dmu4 dmuGBLUP
awk '{print $5,$8}' dmuGBLUP.SOL > temp
mv temp dmuSS.SOL
fi
rm GINV*
rm Gmatrix/gmat.dat
rm Gmatrix/invgmat
rm Gmatrix/gmat
