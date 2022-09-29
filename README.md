# IcelandicCattleBreedingSchemes

This program runs a simulation of a small dairy cattle population undergoing genomic selection.

QMSim was used to simulate the data. The parameters for QMsim are given in the file IceSim.prm. PlinkEdit.sh runs a program that converts the data into plink ped and map format.SelectQTL.R is an R script that converts the files into R objects.

IceSim_main.R is the main file, which is run via the bash file ocs.sh. Input parameters are number of replicate (of the simulated data), method for constructing a relationship matrix, and number of bulls selected in each generation.
Other files are needed for the program to work.
