import sys,os
import os.path
import string
import linecache
import math
import re
import xml.etree.ElementTree as ET

# Currently only capable of creating library for decay calculations.
# No neutron reactions are included in this version.

path = '../ENDF7.1/SubLib_Decay/decay' # set path to point to ENDF7.1 files
outfile1 = open('DeplSubLib_Decay.txt','w')
mass_n = 1.00866491578 # mass of nuetron in amu. Source - ENDF7.1 manual

def Get_ZAID(path,filename):
    file1 = open(path+'/'+filename,'r')
    search = True
    lineNum = 1
    while search:
        line = file1.readline()
        # print(line)
        if lineNum == 2:
            a = line.split()
            ZAtmp = a[0]
            AWRtmp = a[1]
            ZA = re.split("([\+\-])", ZAtmp)
            ZA.insert(1, "E")
            ZA = "".join(ZA[0:])
            # AWR = re.split("([\+\-])", AWRtmp)
            # AWR.insert(1, "E")
            # AWR = "".join(AWR[0:])
            search = False
        lineNum += 1
    # A = float(AWR) * mass_n # use to get Atomic Mass
    # print("A = "+str(A))
    # Z = math.ceil((float(ZA) - A)/1000.) # calculates the atomic number
    # print("Z = "+str(Z))
    return (float(ZA)*10)

for filename in os.listdir(path):
    ZAID = Get_ZAID(path,filename)
    outfile1.write('ZAID = ' + str(ZAID) + '\n')
