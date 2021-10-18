# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:45:50 2021
Program to
@author: au589863
This program takes the input:
    number of replicates
    path to current working directory
Output is:
    ped and map files based on the QMSim files
"""
import time
import datetime
import pandas as pd
import numpy as np
import sys

assert len(sys.argv[1:]) == 2
file = str(sys.argv[1]) + "/r_IceSim/"
rep = str(sys.argv[2])
#file = 'C:/Users/au589863/OneDrive - Aarhus universitet/Documents/phdProject4_Simulations/MoBPS/Marker_based_matrices_branch_11_08_2021/r_IceSim'
#rep = 10
outname = str(sys.argv[1]) + "/data_" + str(rep) + ".ped"
if int(rep) < 10:
    file = file + "/PBLUP_mrk_00" + str(rep) + ".txt"
elif int(rep) > 9:
    file = file + "/PBLUP_mrk_0" + str(rep) + ".txt"

start = time.time()
now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d %H:%M:%S")
print("Program started at :", now)
#assert len(sys.argv[1:]) == 1
#file = str(sys.argv[1])
# load data
df=[]
with open(file, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        df.append(i)
pdf=pd.DataFrame(df[1:])
print(pdf.info())
del df
pdf.columns = ['ID', 'Genotype']
rep={'0':'11', '2':'22', '3': '12','4':'21'}
print("replace values")
pdf['Genotype']=pdf['Genotype'].replace(rep, regex = True)
len(pdf.Genotype[0])
print("change string into columns")
r=pdf.Genotype.apply(lambda x: pd.Series(list(x)).astype(np.int8))

print("insert 5 columns")
r.insert(loc = 0, column = "ID", value = pdf.ID)
r.insert(loc = 1, column = "ID2", value = pdf.ID)
r.insert(loc = 2, column = "col3", value = np.zeros(len(pdf)).astype(np.int8))
r.insert(loc = 3, column = "col4", value = np.zeros(len(pdf)).astype(np.int8))
r.insert(loc = 4, column = "col5", value = np.zeros(len(pdf)).astype(np.int8))
r.insert(loc = 5, column = "col6", value = np.zeros(len(pdf)).astype(np.int8))
del pdf
print("save text")
print(r.info())
np.savetxt(outname,r, fmt='%s', delimiter=" ")
end = time.time()
print("Computations completed in", end - start, "seconds")
now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d %H:%M:%S")
print("Program terminates at :", now)
