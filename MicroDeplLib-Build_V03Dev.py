import sys,os
import os.path
import string
import linecache
import math
import re
import numpy
#packages needed for xml creation and parsing
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

# Currently only capable of creating library for decay calculations.
# No neutron reactions are included in this version.

path = '../ENDF7.1/SubLib_Decay/decay' # set path to point to ENDF7.1 files
outfile1 = open('DeplSubLib_Decay.txt','w')
mass_n = 1.00866491578 # mass of nuetron in amu. Source - ENDF7.1 manual

## Decay mode interpreter - direct from ENDF7.1 manual.
DecMod = ['gamma-ray','Beta-minus','EC/Beta-plus','Isomeric-Trans','Alpha','neutron-emission','Spont-fission','Proton-emission','N/A','N/A','unknown'] # entries 8 and 9 are left out in manual

def Get_Info(current_file,build_flag): #get ZAID and isotope ID
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
    Z = round((float(ZA) - A)/1000.) # calculates the atomic number
    # print("Z = "+str(Z))
    N = round(A-Z) # calculates the number of neutrons
    # print("N = "+str(N))
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

    if build_flag == False: # need to build IsotopeID.txt
        IsotopeID.write(str(ID)+'\t\t'+str(Z)+' '+str(N)+' '+str(A)+'\n')

    return (ZAID,ID,A,Z,N)

def Get_DecayInfo(count, current_file): #get Half_Life, Decay mode(s), Q value(s), and branching ratio(s)
    file1 = open(current_file,'r')
    lineNum = 0
    check = True
    nu=None
    while check:
        lineNum +=1
        line = file1.readline() #read through each line in file starting at line 1
        if " 1452" in line: #signifies start of total number of neutrons from fission - for decay sublibrary, no energy dependence, no prompt/decay. just total
            nuLibStart = lineNum
            tmp = linecache.getline(current_file,(nuLibStart))
            a = tmp.split()
            LNU = a[3] # indicates either tabulated or polynomial representation. for decay lib - SF is always tabulated
            if LNU ==1:
                print('it''s a polynomial? really? ...quitting.')
                sys.exit()
            elif LNU ==2:
                tmp = linecache.getline(current_file,(nuLibStart+2))
                a = tmp.split()
                nu = a[0]

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
                Half_Life=ScientificNotation(a[0])
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
                    Mode[i] = ScientificNotation(ModeTmp)
                    # print(Mode[i])

                    QTmp = a[2]
                    Q[i] = ScientificNotation(QTmp)
                    # print(QTmp)

                    BRTmp = a[4]
                    BR[i] = ScientificNotation(BRTmp)

                    # print(Mode)
                    # print(Q)
                    # print(BR)

                    i += 1
            check = False


    return(Half_Life, Mode, Q, BR, nu)

def ScientificNotation(tmp): # convert ENDF "1.00000+2" to "1.00000E+2"
    if (('+' in tmp) or ('-' in tmp)):
        tmp = re.split("([\+\-])", tmp) #splits ModeTmp into 3 components
        tmp.insert(1, "E") #inserts the E needed for scientific notation
        tmp = "".join(tmp[0:]) #rejoins array of strings into single string
    else:
        pass
    return(tmp)

def TranslateDecayMode(Mode): #translate RTYP numbers to readable decay types
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
            if ((elem-int(elem))>0.): # checks to see if there are multiple decays, e.g. 1.534
                tmp1 = str(elem).split('.') # splits decimal number into two strings, e.g. ['1', '534']
                # print(tmp1)
                tmp2 = list(map(str, tmp1[1])) # splits up '534' into independent strings - e.g. ['5', '3', '4']
                # print(tmp2)
                tmp1[1]=tmp2 # sets tmp2 into tmp1[1] --> tmp1 = ['1', ['5', '3', '4']]
                # print(tmp1)
                for idx,elem1 in enumerate(tmp1):
                    # print('elem1 = '+str(elem1)+' has count of = '+str(elem1.count(elem1)))
                    if elem1.count(elem1)==0: # if elem1 is the ['5', '3', '4'] term
                        dum1 = numpy.chararray(len(elem1),16,True) # initialize dum1
                        for idx2,elem2 in enumerate(elem1): # for each term, '5','3','4'
                            # print(DecMod[int(elem2)])
                            dum1[idx2] = DecMod[int(elem2)] # ID the mod
                    else:
                        dum0 = DecMod[int(elem)]
                        # print(dum0)
                dum = str(dum0)+str(dum1)
            else:
                dum = DecMod[int(elem)]
            TranslatedMode.append(dum)
    return(NumOfModes, TranslatedMode)

def ID_DaughterProducts(Mode,Z,N,SFyields,current_file):
    if Mode == None:
        Daughters=None
        pass
    else:
        Daughters = []#numpy.chararray(len(Mode),6,True)
        for idx,elem in enumerate(Mode):
            if ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-028_Ni_048.endf') and (elem == 2.0)):
                print('  No data for Ni-48 beta+ decay. Need data for Co-48. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-028_Ni_048.endf') and (elem == 7.7)):
                print('  No data for proton-proton decay. Need data for Co-47. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-098_Cf_239.endf') and (elem == 2.0)):
                print('  No data for Bk-239 from Cf-239 beta+ decay. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-098_Cf_256.endf') and (elem == 4.0)):
                print('  No data for Cf-256 alpha decay (and it\'s probability is 1E-9). Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-099_Es_240.endf') and (elem == 4.0)):
                print('  No data for Es-240 alpha decay. Need data for Bk-236. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-099_Es_243.endf') and (elem == 4.0)):
                print('  No data for Es-243 alpha decay. Need data for Bk-239. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-099_Es_258.endf') and (elem == 2.0)):
                print('  No data for Es-243 Beta+ decay. Need data for Cf-258 (which doesn\'t exist?). Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-104_Rf_253.endf') and (elem == 4.0)):
                print('  No data for Rf-253 alpha decay. Need data for No-249. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == '../ENDF7.1/SubLib_Decay/decay/dec-110_Ds_279m1.endf') and (elem == 4.0)):
                print('  No data for Ds-279m1 alpha decay. Need data for Hs-275. Passing...')
                tmp = 'N/A'
                continue
            # print(elem)
            Ztmp = Z; Ntmp = N
            if ((elem-int(elem))>0.): # checks to see if there are multiple decays, e.g. 1.534
                tmp1 = str(elem).split('.') # splits decimal number into two strings, e.g. ['1', '534']
                # print(tmp1)
                tmp2 = list(map(str, tmp1[1])) # splits up '534' into independent strings - e.g. ['5', '3', '4']
                # print(tmp2)
                tmp1[1]=tmp2 # sets tmp2 into tmp1[1] --> tmp1 = ['1', ['5', '3', '4']]
                # print(tmp1)
                for idx,elem1 in enumerate(tmp1):
                    # print('elem1 = '+str(elem1)+' has count of = '+str(elem1.count(elem1)))
                    if elem1.count(elem1)==0: # if elem1 is the ['5', '3', '4'] term
                        dum1 = numpy.chararray(len(elem1),16,True) # initialize dum1
                        for idx2,elem2 in enumerate(elem1): # for each term, '5','3','4'
                            # print(DecMod[int(elem2)])
                            dum1[idx2],Ztmp,Ntmp = Get_Progeny(int(elem2),Ztmp,Ntmp,current_file) # all dumX variables have been translated to their progeny
                            # print('dum1[idx] = '+str(dum1[idx2]))
                    else:
                        dum0,Ztmp,Ntmp = Get_Progeny(int(elem),Ztmp,Ntmp,current_file)
                        # print(dum0)
                dum = str(dum0)+str(dum1)
                    # print(dum)
            else:
                dum,Ztmp,Ntmp = Get_Progeny(int(elem),Ztmp,Ntmp,current_file)
            Daughters.append(dum)
    return(Daughters)

def Get_Progeny(elem,Z,N,current_file):
    IsotopeList = open('IsotopeID.txt','r')
    if elem == 0: # gamma
        pass
    elif elem == 1: # \beta^-
        N-=1; Z+=1
    elif elem == 2: # e.c. / \beta^+
        N+=1; Z-=1
    elif elem == 3: # I.T.
        pass
    elif elem == 4: # alpha
        N-=2; Z-=2
    elif elem == 5: #neutron emission
        N-=1
    elif elem == 6: #SF (NEED TO DEVELOP)
        if current_file in SFyields:
            print(str(current_file)+' do stuff here\n')
        else:
            print('no SFyield info given for '+str(current_file)+'\n')
    elif elem == 7: # Proton-emission
        Z-=1
    elif elem == 10: # unknown
        pass
    else:
        print('I don''t know what decay type I am trying to translate... Quitting.')
        sys.exit()
    # print(str(Z)+' '+str(N))
    for line in IsotopeList:
        a = line.split()
        # print(a)
        if ((str(Z)==a[1]) and (str(N)==a[2])):
            # print('  '+a[0])
            tmp=a[0]
            break # once the daughter is found, break out or else it will cycle through the entire file every single time is searches
        elif (((Z==3) and (N==5)) or ((Z==3) and (N==6))): #if Li-8 or Li-9 (wrong file format, ticket submitted)
            tmp = 'N/A'

    return(tmp,Z,N)


################################################################################
##              START MAIN                                                    ##
################################################################################

## Check if Isotope List needs to be built.
if not os.path.isfile('IsotopeID.txt'): # if the isotopeID list does not exist, create it
    build_flag = False
    IsotopeID = open('IsotopeID.txt','w')
    IsotopeID.write('ID \t\t Z, N, A \n')
    print('Need to build Isotope List. Once finished, please re-exeute script and decay library will be built.')
else:
    build_flag = True

SFyields = os.listdir('../ENDF7.1/SubLib_SFyields/sfy') #gets list of files with Spontaneous fission yields

# Set up xml output information
XML_out = 'DecayData.xml'
file1 = open(XML_out, 'w', encoding = 'utf-8')
root = ET.Element("Decay_Library", Generator = "INL", Name = "General ENDF7.1", Ver = "1.0")
Lib = ET.ElementTree(root)

## Start looping through endf files.
count = 0 # start counter for number of files program runs through
for filename in os.listdir(path):
    count +=1
    # if (filename == 'dec-002_He_003.endf'):
    #     break

# filename = 'dec-005_B_012.endf'; count = 0

    if ((filename == 'dec-003_Li_008.endf') or (filename == 'dec-003_Li_009.endf')):
        pass
    else:
        current_file = path+'/'+filename

    if build_flag == False:
         ZAID,ID,A,Z,N = Get_Info(current_file,build_flag)
         pass
    elif build_flag == True:
        print(filename)
        ZAID,ID,A,Z,N = Get_Info(current_file,build_flag)
        Half_Life, Mode, Q, BR, nu = Get_DecayInfo(count,current_file)
        # print(Mode)
        NumOfModes, TranslatedMode = TranslateDecayMode(Mode)
        Daughters = ID_DaughterProducts(Mode,Z,N,SFyields,current_file)
        # print('Daughters are: '+str(Daughters)+'\n')
        # outfile1.write('ID = ' + str(ID) + '\t\t\tZAID = ' + str(ZAID) + '\t\tT_{1/2} = ' + str(Half_Life) + '\t\tNum of Decay Modes = ' + str(NumOfModes) + '\n')
        # outfile1.write('\t\t Decay Mode(s) = ' + str(TranslatedMode) + '\n\t\t Q-value(s) = ' + str(Q) + '\n\t\t Branching Ratio(s) = ' + str(BR) + '\n'+'\t\t Daughter Product(s) = '+str(Daughters)+'\n\n')

        isotope = ET.SubElement(root, ("Isotope"), Halflife = str(Half_Life), Name = str(ID), ZAID = str(ZAID) )
        Type = ET.SubElement(isotope, "type")
        Type.text = str(TranslatedMode).strip("[]")
        Progeny = ET.SubElement(isotope, "daughters")
        Progeny.text = str(Daughters).strip("[]")
        Branch_Ratio = ET.SubElement(isotope, "branch_ratio")
        Branch_Ratio.text = str(BR).strip("[]")


Lib.write(file1, encoding='unicode')
file1.close

### temp workaround for parsing....
# os.system('xmllint --format DecayData.xml > DecayData_tmp.xml')
# os.system('rm DecayData.xml')
# os.system('mv DecayData_tmp.xml DecayData.xml')

# doc2 = minidom.parse(XML_out)

#### The code below works on simple toy problems but not here.... so use workaround system commands above.
## use xml.dom.minidom and parse out previously created xml file.
# xml = minidom.parse('DecayData.xml')
# pretty_xml_as_string = xml.toprettyxml()
# file2 = open('DecayData.xml', 'w', encoding = 'utf-8')
# file2.write(pretty_xml_as_string)
# file2.close()
