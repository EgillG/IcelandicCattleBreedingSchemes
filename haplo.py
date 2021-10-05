# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:58:37 2020
try to fix for loop
@author: au589863
"""
import os
import sys
import datetime
import time
import numpy as np
import pandas as pd

start = time.time()
now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d %H:%M:%S")
print("Program started at :", now)
assert len(sys.argv[1:]) == 1
method = str(sys.argv[1])
# load data
file="haplo.txt"
df=[]
with open(file, 'r') as fp:
     for cnt, line in enumerate(fp):
         i=line.rstrip("\n").split("\t")
         df.append(i)
# make pandas dataframe, set first row as columns
cols=df[0]
pddf=pd.DataFrame(df[1:])
del df
pddf.columns=cols
# info about dataframe
pddf.info()
# combine all columns into a string.
haplo = pddf.agg(''.join, axis=0)
# Make dataframe out of the genotype string:
hapdf=pd.DataFrame(haplo, columns = ["haplo1"])
del haplo
# make a column of IDs, with the chromosome number removed (last 5 characters)
hapdf['ID_chr'] = hapdf.index.astype(str)
hapdf['ID_chr'] = hapdf['ID_chr'].str.slice(0,-5)
hapdf

# make a column of chromosome number, haplo1 and haplo2
hapdf['sort'] = pd.Series(['haplo' + str(i) for i in [1,2]*int(len(hapdf.haplo1)/2)]).values
# use pivot attribute to turn rows into columns, returning
# a dataframe with one row for each animal, and one column for each haplotype.
hapdf = hapdf.pivot(columns = 'sort', values='haplo1',  index='ID_chr')
#Haplotype lengths to consider. Should be an input to the program.
j=10

# make list for storing dataframes
lst=[]
# Gera dictionary sem tekur inn haplótýpurnar og býr til raðnúmer í staðinn.
# Síðan er hægt að raða dálkunum samkvæmt raðnúmerinu
# Ef raðnúmerin eru gerð innan hverrar haplótýpubreiddar er sorteringin ekki 
# tímafrek.

# Loop over all the markers in steps of j.
numAn=len(hapdf)
numMarkers=len(hapdf.haplo1[0])
# make a dataframe with all haplotyes, for allele coding of the pseudo SNP.
allhap=pd.concat([hapdf.haplo1,hapdf.haplo2])
for i in range(0,numMarkers,j):
    # for each j markers, use cat.codes to recode the unique haplotypes.
    # Use cat.codes to assign the haplotype string (01010)
    # to a pseudo SNP name in a dictionary.
    haplos = pd.concat([allhap.str.slice(i,i+j).astype(str),(allhap.str.slice(i,i+j).astype('category').cat.codes.astype(str))], axis=1)
    haplos.columns=['hap','code']
    # Make a dictionary using list comprehension
    dictionary = dict([(a,int(b)) for a,b in zip(haplos.hap,haplos.code)])
    # replace genotype codes with new integer code generated with cat.codes above.
    haplo1 = hapdf['haplo1'].str.slice(i,i+j).astype(str).replace(dictionary) 
    haplo2 = hapdf['haplo2'].str.slice(i,i+j).astype(str).replace(dictionary)
    count=1
    # Here I have to loop through all the haplotypes, and make a list of positions
    # which will be used to insert columns of zeros if the allele is not present at
    # all in either haplo1 or hapl2. 
    # zer is the list of index positions where a column of zeros is inserted later.
    zer = []
    # dictionary with string name as key and integer number as value.
    # This is so that it is possible to order the columns for input to GMATRIX,
    # where two consecutive columns represent parental chromosome, and 
    # each row has 1 or 0 indicating whether the allele is present or not in 
    # the animal.
    indi={}
    # count variable is to assign integer column names
    count=1
    # f is a list of sorted dictionary values, which is iterated over
    f=list(dictionary.values())
    f.sort()
    for bla in f:
        if bla in haplo1.unique() and bla in haplo2.unique():
#            print("in haplo1 and in haplo2")
            indi['haplo1_' + str(bla) + '_1'] = count
            indi['haplo2_' + str(bla) + '_2'] = count + 1
            count = count + 2
        elif bla in haplo1.unique() and bla not in haplo2.unique():
#            print("in haplo1 but not in haplo2")
            indi['haplo1_' + str(bla) + '_1'] = count
            indi['haplo2_' + str(bla) + '_2'] = count + 1
            zer.append(count + 1)
            count = count + 2
        elif bla in haplo2.unique() and bla not in haplo1.unique():
#            print("in haplo2 but not in haplo1")
            indi['haplo1_' + str(bla) + '_1'] = count
            indi['haplo2_' + str(bla) + '_2'] = count + 1
            zer.append(count)
            count = count + 2
    # make new variables haplo1 and haplo2. They are used to 
    haplo1 = hapdf['haplo1'].str.slice(i,i+j).astype(str).replace(dictionary).astype(str) + '_1'
    haplo2 = hapdf['haplo2'].str.slice(i,i+j).astype(str).replace(dictionary).astype(str) + '_2'
    # make list for column of zeros
    zeros = [0]*len(haplo1)
    # concatenate haplo1 and haplo2 along columns (axis = 1)
    # Then I have a dataframe with ID of animal in rows, and parental haplotype alleles
    # as columns for the window considered.
    # Values are haplotype allele code, underscore, 1 or 2 designate chromosome.
    conc=pd.concat([haplo1, haplo2], axis=1)
    # use get_dummies to get dummy variables for these pseudo SNPs.
    dummy=pd.get_dummies(conc)
    # dummy.sum().sum() gives the number of animals in the relationship matrix.
    # This short loop inserts columns of zeros where either chromosome (set) 1, or 2, 
    # was not present above
    for idx in zer:
        dummy.insert(loc = idx-1, column = idx, value = zeros)
    # rename the columns according to the indi dictionary, allowing sort.
    dummy = dummy.rename(columns = indi)
    # sort the columns
    dummy=dummy.sort_index(axis=1)
    dropList=[]
#    for i in range(1,len(dummy.columns+1), 2):
#        if dummy.loc[:,[i,i+1]].sum().sum()/numAn <0.1:
#            dropList.append(i)
#            dropList.append(i+1)
#    dummy.drop(dropList, axis = 1,inplace = True)    
    # append to lst
    lst.append(dummy.astype(np.int8))#lst.append(dummy.rename(columns = {value : key for (key, value) in indi.items()}))
# At last, combine the dataframes in one big dataframe using pandas concat
dfOut=pd.concat(lst, sort=True,axis=1) #ignore_index = True?
dfOut=dfOut+1
# Write a file with IDs and new IDs
#IDout = pd.concat([pd.Series(np.arange(1, len(dfOut) + 1)),pd.Series(dfOut.index),pd.Series(np.ones([len(dfOut)])).astype(np.int8)], axis = 1)
IDout = pd.concat([pd.Series(np.arange(1, len(dfOut) + 1)),pd.Series(np.arange(1, len(dfOut) + 1)),pd.Series(np.ones([len(dfOut)])).astype(np.int8)], axis = 1)

IDout.to_csv("gmatrix.id", sep = " ", index = False, header = False)
numPsSnp = len(dfOut.columns)
mapfile = pd.concat([pd.Series(np.ones(numPsSnp)).astype(np.int8),pd.Series(np.ones(numPsSnp)).astype(np.int8),pd.Series(np.ones(numPsSnp)).astype(np.int8)], axis = 1)
mapfile.to_csv("gmatrix.map", sep = " ", index = False, header= False)
   
# Write data
dfOut=dfOut.astype(str).agg(''.join, axis=1)
dfOut=pd.DataFrame(dfOut)
dfOut.insert(loc = 0, column = "ID", value = np.arange(1, len(dfOut) + 1))
Savename='gmatrix.dat'
np.savetxt(Savename,
            dfOut,fmt='%s', delimiter=" ")

end = time.time()
print("Computations completed in", end - start, "seconds")
now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d %H:%M:%S")
print("Program terminates at :", now)
