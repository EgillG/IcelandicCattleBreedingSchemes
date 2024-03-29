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
[egill@has802 QMSimTest]$ cat QMsimPar.py
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 12:20:11 2021

@author: au589863
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:44:25 2021

@author: au589863
"""
import sys
import os
assert len(sys.argv[1:]) == 5
NumRep = str(sys.argv[1]) # Number of replicates
NumChrom = str(sys.argv[2]) # Number of chromosomes
NumMark = str(sys.argv[3]) # Number of markers on each chromosome
SizeChrom = str(sys.argv[4]) # Size of each chromosome (centiMorgans)
wd = str(sys.argv[5]) # Current working drive
#wd =  'C:\\Users\\au589863\OneDrive - Aarhus universitet\Documents\phdProject4_Simulations\MoBPS\8Test_22_06_2021'
#NumRep = 10
#NumChrom = 29
#NumMark = 1700
#SizeChrom = 150
#heritability = 0.4
par = '/* \n\
Builds historical population\n\
\n\
*/\n\
\n\
/*******************************\n\
 **     Global parameters     **\n\
 *******************************/\n\
title = "Simulation of Icelandic Dairy Cattle";\n\
nrep  =' + str(NumRep) + ';                    //Number of replicates\n\
h2    = 0.4;                  //Heritability\n\
qtlh2 = 0.4;                  //QTL heritability\n\
phvar = 1.0;                  //Phenotypic variance\n\
nthread = 2;\n\
no_male_rec;\n\
/*******************************\n\
 **   Historical population   **\n\
 *******************************/\n\
// Use the same as in other papers\n\
begin_hp;\n\
   hg_size = 1000 [0]          //Size of the historical generations\n\
             1000 [3000];\n\
   nmlhg   = 500;              //Number of males in the last generation\n\
end_hp;\n\
\n\
/*******************************\n\
 **        Populations        **\n\
 *******************************/\n\
\n\
begin_pop = "PBLUP";\n\
   begin_founder;\n\
      male   [n = 500, pop = "hp"];\n\
      female [n = 500, pop = "hp"];\n\
   end_founder;\n\
   ls  = 1;                 //Litter size\n\
   ng  = 1;                 //Number of generations\n\
   md  = rnd;               //Mating design\n\
   sr  = 1.0;               //Replacement ratio for sires: Discrete generations\n\
   dr  = 0.5;               //Replacement ratio for dams:  50% replaced\n\
   cd  = age;\n\
   pmp = 0.5 /fix;          //Proportion of male progeny\n\
   begin_popoutput;\n\
        data;\n\
        stat;\n\
        genotype /snp_code /gen 0;\n\
        allele_freq;\n\
        ld;\n\
   end_popoutput;\n\
end_pop;\n\
\n\
/*******************************\n\
 **          Genome           **\n\
 *******************************/\n\
begin_genome;'
genome = ""
for i in range(1,int(NumChrom)+1):
    genome = genome + ' begin_chr = ' + str(1) + ';                                                     //Chromosome id\n\
        chrlen = ' + str(SizeChrom) + ';                                                //Chromosome lenght\n\
        nmloci =' + str(NumMark) + ';                                           //Number of marker loci on the chromosome\n\
        mpos   = even;                                          //Markers positions\n\
        nma    = all 2;                                         //Number of marker alleles\n\
        maf    = eql;                                           //Initial frequency of marker alleles\n\
        nqloci = 500;                                           //Number of QTL loci\n\
        qpos   = rnd;                                           //QTLs positions\n\
        nqa    = all 2;                                         //Number of QTL alleles\n\
        qaf    = eql;                                           //Initial frequency of QTL alleles\n\
        qae    = rndg 0.4;                                      //Distribution of QTL allele effects\n\
     end_chr;'
last = '    mmutr     = 2.5e-5 /recurrent;      //Marker mutation rate\n\
   qmutr     = 0;                                       //QTL mutation rate\n\
end_genome;\n\
\n\
\n\
/*******************************\n\
 **       Output options      **\n\
 *******************************/\n\
begin_output;\n\
   linkage_map;\n\
end_output;\n\
\n\
/ \n\
'

out=open(wd+"/IceSim.prm", 'w')
L = [par+genome+last]
out.writelines(L)
out.close()
