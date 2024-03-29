# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:58:24 2021

This program reads a file called SelCands, which contains the ID, sexc
and number of progeny of animals.
If number of progeny is higher than 0,
then the ID and sex is stored in a dictionary.

Then the program reada relationships from the file
ReductedMatrix. Male-male relationshpis are written to a file
EVArels.txt, but female to male relationships are added and stored.
These sums are divided with number of female selection candidates to
obtain the average relationship of females to males.

Lastly, all pairwise female-female relationships are summed and divided
with number of relationships to obtain the average relationship among
females.
The pseudo-female is then constructed and added to EVArels.txt
@author: au589863
"""
import datetime
import time

start = time.time()
now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d %H:%M:%S")
print("G_matrixPreparation_3 program started at :", now)

print("Reading SelCands file")
DicFile = "SelCands"
dictionary = {}
with open(DicFile, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        dictionary[i[0]] = int(i[1])

MaleList = {i:0 for i in dictionary if dictionary[i]==1}
numFemales = len({i:0 for i in dictionary if dictionary[i]==2})
numMales = len(MaleList)
file = "ReducedMatrix"
df=[]
print("Number of males in SelCands is: " + str(numMales) + ". Number of females is: " + str(numFemales) +".")

print("Reading ReducedMatrix and writing relationships to EVArels.txt")
count=0
NewFile = open("EVArels.txt","w")
# Read male -female relationships
with open(file, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        if i[0] not in dictionary or i[1] not in dictionary:
            continue
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
print("Read and wrote " + str(count) + " male-to-male relationships to EVArels.txt.")
meanRel = 0
count=0
with open(file, 'r') as fp:
    for cnt, line in enumerate(fp):
        i=line.rstrip("\n").split( )
        if i[0] not in dictionary or i[1] not in dictionary:
            continue
        elif dictionary[i[0]] == 2 and dictionary[i[1]] == 2:
                meanRel=meanRel + float(i[2])
                count +=1 
print("Read " + str(count) + " female-to-female relationships.")
print("Average relationships of females to females was: " + str(meanRel/count) +".")

assert count == (numFemales*(numFemales-1))/2 + numFemales
meanRel/count
# Testing code
# sum([float(i[2]) for i in df])/count
NewFile = open("EVArels.txt","a")
print("Writing female relationships to EVArels.txt")
for line in List:
    NewFile.write(str(line[0]) + "\t00001\t" + str(line[1]) +"\n")
NewFile.write("00001\t00001\t"+str(meanRel/count)+ "\n")
NewFile.close()

end = time.time()
print("Computation of relationships completed in", end - start, "seconds")
now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d %H:%M:%S")
print("G_matrixPreparation_3 program terminates at :", now)
