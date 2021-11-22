method=${1?Error: method not given to DMU.sh}
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

cd $SLURM_SUBMIT_DIR/${method}/
rm -rf /scratch/$USER/$SLURM_JOBID

echo "Gmatrix Job completed at $(date '+%y-%m-%d %H:%M:%S')"
cut -f 1 -d ' ' Gmatrix/gmat.dat > genotyped
echo running single step GBLUP DMU
bash r_dmu5 dmuSS
awk '{print $5,$8}' dmuSS.SOL > temp
mv temp dmuSS.SOL
