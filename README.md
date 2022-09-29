# IcelandicCattleBreedingSchemes

This program runs a simulation of a small dairy cattle population undergoing genomic selection.

QMSim was used to simulate the data. The parameters for QMsim are given in the file IceSim.prm. PlinkEdit.sh runs a program that converts the data into plink ped and map format.SelectQTL.R is an R script that converts the files into R objects.

IceSim_main.R is the main file, which is run via the bash file ocs.sh.
EvaPar.py and gmatpar.py are programs that write parameter files for the programs EVA and GMATRIX.
EVA.sh is a bash script that takes arguments from the main file.
