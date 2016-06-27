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

def Get_Info(current_file):
    ################################################################
    ##           get base isotope ZAID                            ##
    ################################################################
    tmp = linecache.getline(current_file,2,None) #pulls line 2 from current_file into cache
    line2 = tmp.split() #splits line into array of strings
    ZAtmp = line2[0] #gets ZA from line 2
    ZA = re.split("([\+\-])", ZAtmp) #splits ZA into 3 components
    ZA.insert(1, "E") #inserts the E needed for scientific notation
    ZA = "".join(ZA[0:]) #rejoins array of strings into single string
    ZAID = int(float(ZA)*10.0)

    ################################################################
    ####    get Isotope Name                                      ##
    ################################################################
    tmp = linecache.getline(current_file,6,None)
    line6 = tmp.split()
    IDtmp = line6[0]
    IDtmp = re.split("-",IDtmp)
    IDtmp = IDtmp[1]
    ## the following code block can be used to explicitely calculate Z and A
    AWRtmp = line2[1]
    AWR = re.split("([\+\-])", AWRtmp)
    AWR.insert(1, "E")
    AWR = "".join(AWR[0:])
    A = float(AWR) * mass_n # use to get Atomic Mass
    # print("A = "+str(round(A)))
    # Z = math.ceil((float(ZA) - A)/1000.) # calculates the atomic number
    # print("Z = "+str(Z))
    ID = IDtmp+str(round(A))

    ################################################################
    ##      check for metastable state                            ##
    ################################################################
    tmp = linecache.getline(current_file,3,None) #pulls line 3 into cache
    line3 = tmp.split() #split line into array of strings
    if int(line3[3]) > 0:
        ZAID = str(int(ZAID/10))+str('1') #'1' is used for metastable state
        ID = ID+str('m')

    ################################################################
    ####    get half life (seconds)                               ##
    ################################################################
    

    linecache.clearcache()
    return (ZAID,ID)

for filename in os.listdir(path):
    current_file = path+'/'+filename
    ZAID,ID = Get_Info(current_file)
    outfile1.write('ID = ' + str(ID) + '\t\t\tZAID = ' + str(ZAID) + '\n')
