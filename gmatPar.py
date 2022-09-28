# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 10:21:07 2020

@author: au589863
"""
import sys
import os
assert len(sys.argv[1:]) == 3
pwd = str(sys.argv[1]) + "/Gmatrix"
method = str(sys.argv[2])
singleS = str(sys.argv[3])
print(pwd)
outname= pwd + "/gmat.par"
if not os.path.exists("Gmatrix"):
    os.mkdir("Gmatrix")
if method[1] == "1":
    scalemtd = 1
elif method[1] == "2":
    scalemtd = 2
    
if len(method) > 2:
    freqmtd = 11
else:
    freqmtd = 1
if singleS == "Yes":
   out_invG = str(0)
   out_G = str(1)
elif singleS == "No":
   out_invG = str(1)
   out_G = str(0)

par = "$MAPFILE \n" + '\"' + pwd + "/gmat.map\" \n\
$MARKERFILE \n" + '\"' + pwd + "/gmat.dat\" \n\
$MAPBASEFILE \n" + '\"' + pwd + "/mapbase.dat\" \n\
$IDFILE \n" + '\"' + pwd + "/gmat.id\" \n\
$GMATRIXFILE \n" + '\"' + pwd + "/gmat\" \n\
$IGMATRIXFILE \n" + '\"' + pwd + "/invgmat\"\n\n\
%PARSTART\n\n\
$DATATYPE\n\
2\n\
$DELIMITER\n\
1\n\
$NSKIP\n\
1\n\
$A1A1CODE\n\
0\n\
$MISSGENOTYPE\n\
9\n\
$MINMAF\n\
0.00000001\n\
$FREQMETHOD\n" + \
str(freqmtd) + "\n\
$SCALEMETHOD\n" + \
str(scalemtd) + "\n\
$DIAG_ADD\n\
0.01\n\
$G_ADD\n\
0.00\n\
$DIAG_ONE\n\
0\n\
$AGSCALE\n\
0\n\
$PROP_A_to_G\n\
0\n\
$ADJSCALEMISS\n\
0\n\
$CAL_DET\n\
1\n\
$OUTPRECISION\n\
4\n\
$OUT_GMATRIX\n" + \
out_G + "\n\
$OUT_IGMATRIX\n" +\
out_invG + "\n\
$OUT_AMATRIX\n\
0\n\
$OUT_AMATRIX11\n\
0"

out=open(outname, 'w')
L = [par]
out.writelines(L)
out.close()
