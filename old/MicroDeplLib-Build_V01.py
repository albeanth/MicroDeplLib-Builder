import sys,os
import os.path
import string
import linecache
import math
import re
import numpy
import xml.etree.ElementTree as ET

# Currently only capable of creating library for decay calculations.
# No neutron reactions are included in this version.

path = '../ENDF7.1/SubLib_Decay/decay' # set path to point to ENDF7.1 files
outfile1 = open('DeplSubLib_Decay.txt','w')
mass_n = 1.00866491578 # mass of nuetron in amu. Source - ENDF7.1 manual

## Decay mode interpreter - direct from ENDF7.1 manual.
DecMod = ['gamma-ray','Beta-minus','EC/Beta-plus','Isomeric-Trans','Alpha','neutron-emission','Spont-fission','Proton-emission','N/A','N/A','unknown'] # entries 8 and 9 are left out in manual

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

    linecache.clearcache()
    return (ZAID,ID)

def Get_DecayInfo(count, current_file):
    file1 = open(current_file,'r')
    lineNum = 0
    check = True
    while check:
        lineNum +=1
        line = file1.readline() #read through each line in file starting at line 1
        if " 8457" in line: #signifies start of library 8 (decay) MT=457, rad decay data set
            DecLibStart = lineNum
            # get half life
            tmp = linecache.getline(current_file,(DecLibStart+1),None)
            a = tmp.split()
            Half_Life=a[0]
            if Half_Life == '0.000000+0': #signifies that isotope is stable no decay data
                Half_Life = str('Infinity')
                Mode = None
                Q = None
                BR = None
                break
            else:
                # get number of decay modes (NDK)
                tmp = linecache.getline(current_file,(DecLibStart+3),None)
                a = tmp.split()
                # print(a)
                NDK=a[5]
                # print(NDK)
                if ((count>999) & (NDK.endswith(str(count)))):
                    NDK = int(NDK[:-len(str(count))])
                else:
                    NDK=int(NDK)
                # get decay info: Modes, Q values, and Branching Ratios (BR)
                Mode = numpy.empty(NDK,dtype=float) #create decay mode vector of length equal to NDK
                Q = numpy.empty(NDK,dtype=float) #create Q value vector for each mode
                BR = numpy.empty(NDK,dtype=float) #create Branching ratio  vector of length equal to NDK
                RTYP_start = DecLibStart+4 # line number in which mode of decay identifiers (RYTP) start
                i = 0
                while (i<NDK):
                    tmp = linecache.getline(current_file,(RTYP_start+i),None)
                    a = tmp.split()
                    # print(a)
                    ModeTmp = a[0]
                    ModeTmp = re.split("([\+\-])", ModeTmp) #splits ModeTmp into 3 components
                    ModeTmp.insert(1, "E") #inserts the E needed for scientific notation
                    ModeTmp = "".join(ModeTmp[0:]) #rejoins array of strings into single string
                    Mode[i] = float(ModeTmp)
                    # print(Mode[i])

                    QTmp = a[2]
                    QTmp = re.split("([\+\-])", QTmp) #splits ModeTmp into 3 components
                    QTmp.insert(1, "E") #inserts the E needed for scientific notation
                    QTmp = "".join(QTmp[0:]) #rejoins array of strings into single string
                    Q[i] = float(QTmp)
                    # print(QTmp)

                    BRTmp = a[4]
                    BRTmp = re.split("([\+\-])", BRTmp) #splits ModeTmp into 3 components
                    BRTmp.insert(1, "E") #inserts the E needed for scientific notation
                    BRTmp = "".join(BRTmp[0:]) #rejoins array of strings into single string
                    BR[i] = float(BRTmp)

                    # print(Mode)
                    # print(Q)
                    # print(BR)

                    i += 1
            check = False
    return(Half_Life, Mode, Q, BR)


################################################################################
##              START MAIN                                                    ##
################################################################################
count = 0 # start counter for number of files program runs through
for filename in os.listdir(path):
    count +=1
    # if (filename == 'dec-002_He_010.endf'):
    #     break
# filename = 'dec-039_Y_091m1.endf'
# count = 1000

    print(filename)

    if ((filename == 'dec-003_Li_008.endf') or (filename == 'dec-003_Li_009.endf')):
        pass
    else:
        current_file = path+'/'+filename
        ZAID,ID = Get_Info(current_file)
        Half_Life, Mode, Q, BR = Get_DecayInfo(count,current_file)
        # print(Mode)
        if Mode == None:
            TranslatedMode = None
            NumOfModes = None
        else:
            NumOfModes=len(Mode)
            TranslatedMode = numpy.chararray((NumOfModes,7),16,True) # 16 is the length of the longest element in DecMod. True - changes unicode setting to true
            # print(Mode)
            # print(TranslatedMode)
            TranslatedMode=[]
            for idy,elem in enumerate(Mode):
                # print("base element = "+str(elem))
                if ((elem-int(elem))>0.):
                    tmp1 = str(elem).split('.')
                    # print(tmp1)
                    tmp2 = list(map(str, tmp1[1]))
                    # print(tmp2)
                    tmp1[1]=tmp2
                    # print(tmp1)
                    for idx,elem1 in enumerate(tmp1):
                        # print('elem1 = '+str(elem1)+' has count of = '+str(elem1.count(elem1)))
                        if elem1.count(elem1)==0:
                            dum1 = numpy.chararray(len(elem1),16,True)
                            for idx2,elem2 in enumerate(elem1):
                                # print(DecMod[int(elem2)])
                                dum1[idx2] = DecMod[int(elem2)]
                        else:
                            dum0 = DecMod[int(elem)]
                            # print(dum0)
                    dum = str(dum0)+str(dum1)
                else:
                    dum = DecMod[int(elem)]
                TranslatedMode.append(dum)
            # print(TranslatedMode)


    outfile1.write('ID = ' + str(ID) + '\t\t\tZAID = ' + str(ZAID) + '\t\tT_{1/2} = ' + str(Half_Life) + '\t\tNum of Decay Modes = ' + str(NumOfModes) + '\n')
    outfile1.write('\t\t Decay Mode(s) = ' + str(TranslatedMode) + '\n\t\t Q-value(s) = ' + str(Q) + '\n\t\t Branching Ratio(s) = ' + str(BR) + '\n\n')
