# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:58:24 2021

Here I start by making a dictionary for selection candidates.
Then the program read female - male relationships only, and adds these values
to the dictionary. After all values have been stored, the numbers
are divided with the number of female selection candidates.

@author: au589863
"""

DicFile = "SelCands"
dictionary = {}
with open(DicFile, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        if int(i[2]) > 0:
            dictionary[i[0]] = int(i[1])

MaleList = {i:0 for i in dictionary if dictionary[i]==1}
numFemales = len(dictionary) - len(MaleList)
numMales = len(MaleList)
file = "ReducedMatrix"
df=[]


count=0
NewFile = open("EVArels.txt","w")
# Read male -female relationships
with open(file, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        if dictionary[i[0]] == 1 and dictionary[i[1]] == 1:
                 NewFile.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
                 count += 1
        if dictionary[i[0]] == 1 and dictionary[i[1]] == 2:
            # i.append(dictionary[i[0]])
            # i.append(dictionary[i[1]])
            # df.append(i)
            MaleList[i[0]] = MaleList[i[0]] + float(i[2])
        elif dictionary[i[0]] == 2 and dictionary[i[1]] == 1:
            # i.append(dictionary[i[0]])
            # i.append(dictionary[i[1]])
            # df.append(i)
            MaleList[i[1]] = MaleList[i[1]] + float(i[2])
assert count == (numMales*(numMales-1))/2+numMales
NewFile.close()
List = [[k,float(v)/numFemales] for k,v in MaleList.items()]
len(List)

# Here is test for the means
# import pandas as pd
# dataF=pd.DataFrame(df)
# dataF.columns = ["A","B","C","D","E"]
# dataF[["C"]] = dataF[["C"]].apply(pd.to_numeric)
# group_df=dataF[["B","C"]].groupby("B").mean()
# List[:5]

Inbreeding = 0
InbreedCount = 0
meanRel = 0
count=0
with open(file, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        if dictionary[i[0]] == 2 and dictionary[i[1]] == 2:
            if i[0] != i[1]:
                meanRel=meanRel + float(i[2])
                count +=1 
#                df.append(i)
            if i[0] == i[1]:
                Inbreeding=Inbreeding+float(i[2])
                InbreedCount += 1
assert InbreedCount == numFemales
assert count == (numFemales*(numFemales-1))/2
Inbreeding/InbreedCount
meanRel/count
# Testing code
# sum([float(i[2]) for i in df])/count
NewFile = open("EVArels.txt","a")
for line in List:
    NewFile.write(str(line[0]) + "\t00001\t" + str(line[1]) +"\n")
NewFile.write("00001\t00001\t"+str(meanRel/count)+ "\n")
NewFile.close()
