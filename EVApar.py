# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:44:25 2021

@author: au589863
"""
import sys
import os
assert len(sys.argv[1:]) == 3
method = str(sys.argv[1])
numBulls = str(sys.argv[2])
wd = str(sys.argv[3])
fileG = 'file'
outname= wd + "/" + "/EVA.prm"
if not os.path.exists("Gmatrix"):
    os.mkdir("Gmatrix")
if not os.path.exists("_evaSim"):
    os.mkdir("_evaSim")

par = "&DataParameters \n\
   DataFile                     = '" + wd + "/evaIn.txt' ,\n\
   ResultsDirectory             = '" + wd + "/evaSim/" + "' ,\n\
   IgnoreParentalPedigreeErrors = .true. , \n\
   RecodeFile                   = 'RecodedIds.dat'\n\
/ \n\
&PopulationHistory \n\
  PCI_ngen    = 5\n\
/ \n\
&RelationshipMatrix \n\
  Source  = " + fileG + ",\n\
  Gfile = '" + wd + "/EVArels.txt',\n\
  TimeSteps = 0\n\
/ \n\
&OCSParameters \n\
  nmatings         = " +numBulls + ",\n\
  optimise         = 'constraint',\n\
  dFconstraint     = 0.005,\n\
  LimitMaleMatings = 1,\n\
/ \n\
&AlgorithmParameters \n\
  generations                   = 10000 ,\n\
  popsize                       = 100 ,\n\
  n_offspring                   = 50 ,\n\
  restart_interval              = 2000 ,\n\
  exchange_algorithm            = 500,\n\
  mutate_probability            = 0.001 ,\n\
  crossover_probability         = 0.20,\n\
  directed_mutation_probability = 0.01,\n\
  seed_rng                      = 0 \n\
/ \n\
\n\
&MatingOptions\n\
  MatingStrategy='random' ,\n\
  RepeatedMatings=.true.\n\
/ "

out=open(outname, 'w')
L = [par]
out.writelines(L)
out.close()
