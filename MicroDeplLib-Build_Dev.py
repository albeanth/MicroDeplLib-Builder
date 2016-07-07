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

import translators

dec_path = '../ENDF7.1/SubLib_Decay/decay' # set path to point to ENDF7.1 files
nRxn_path = '../ENDF7.1/SubLib_NeutronRxn/neutrons' # set path to point to ENDF7.1 files
mass_n = 1.00866491578 # mass of nuetron in amu. Source - ENDF7.1 manual

dec_List = os.listdir(dec_path) # get list of isotope files in decay sublibrary
nRxn_List = os.listdir(nRxn_path) # get list of isotope files in neutron reactions

## GENERAL PURPOSE FUCNTIONS
def ScientificNotation(tmp): # convert ENDF "1.00000+2" to "1.00000E+2"
    if (('+' in tmp) or ('-' in tmp)):
        tmp = re.split("([\+\-])", tmp) #splits ModeTmp into 3 components
        tmp.insert(1, "E") #inserts the E needed for scientific notation
        tmp = "".join(tmp[0:]) #rejoins array of strings into single string
    else:
        pass
    return(tmp)

def Get_Info(current_file): #get ZAID and isotope ID
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

    return (ZAID,ID,Z,N)

def NuclideIdentifier(Z,N):
    IsotopeList = open('IsotopeID.txt','r')
    for line in IsotopeList:
        a = line.split()
        if ((str(Z)==a[1]) and (str(N)==a[2])):
            name=a[0]
            break # once the daughter is found, break out or else it will cycle through the entire file every single time is searches
        elif (((Z==3) and (N==5)) or ((Z==3) and (N==6))): #if Li-8 or Li-9 (wrong file format, ticket submitted)
            name = 'N/A'
    return(name)

## DECAY SPECIFIC FUNCTIONS
def Get_DecayInfo(count, current_file): #get Half_Life, Decay mode(s), Q value(s), and branching ratio(s)
    file1 = open(current_file,'r')
    lineNum = 0
    nu=None # number of neutrons per fission event
    for line in file1:
        lineNum+=1
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
                # print(a)
                nu = a[0]

        elif " 8457" in line: #signifies start of library 8 (decay) MT=457, rad decay data set
            DecLibStart = lineNum
            # print(DecLibStart)
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
                # Mode = [] # this style might clean up code a bit. but it throws errors due to some format issues.
                # Q = []
                # BR = []
                RTYP_start = DecLibStart+4 # line number in which mode of decay identifiers (RYTP) start
                i = 0
                # print(NDK)
                while (i<NDK):
                    tmp = linecache.getline(current_file,(RTYP_start+i),None)
                    a = tmp.split()
                    # print(a)

                    ModeTmp = a[0]
                    Mode[i] = ScientificNotation(ModeTmp)
                    # Mode.append(ScientificNotation(ModeTmp))
                    # print(Mode[i])

                    QTmp = a[2]
                    Q[i] = ScientificNotation(QTmp)
                    # Q.append(ScientificNotation(QTmp))
                    # print(QTmp)

                    BRTmp = a[4]
                    BR[i] = ScientificNotation(BRTmp)
                    # BR.append(ScientificNotation(BRTmp))

                    # print(Mode)
                    # print(Q)
                    # print(BR)

                    i += 1
                break
        # elif '8  099999' in line:
    return(Half_Life, Mode, Q, BR, nu)

def TranslateDecayMode(Mode,Z,N,SFyields,current_file): #translate RTYP numbers to readable decay types
    if Mode == None:
        TranslatedMode = None
        NumOfModes = None
        Daughters = None
    else:
        NumOfModes=len(Mode)
        Daughters = []
        TranslatedMode=[]
        for idy,elem in enumerate(Mode):
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
            Ztmp = Z; Ntmp = N
            if ((elem-int(elem))>0.): # checks to see if there are multiple decays, e.g. 1.534
                tmp1 = str(elem).split('.') # splits decimal number into two strings, e.g. ['1', '534']
                tmp2 = list(map(str, tmp1[1])) # splits up '534' into independent strings - e.g. ['5', '3', '4']
                tmp1[1]=tmp2 # sets tmp2 into tmp1[1] --> tmp1 = ['1', ['5', '3', '4']]
                for idx,elem1 in enumerate(tmp1):
                    if elem1.count(elem1)==0: # if elem1 is the ['5', '3', '4'] term
                        DecNameL = numpy.chararray(len(elem1),16,True) # initialize dum1
                        ProgNameL = numpy.chararray(len(elem1),16,True) # initialize dum1
                        for idx2,elem2 in enumerate(elem1): # for each term, '5','3','4'
                            DecNameL[idx2],Ztmp,Ntmp = translators.Progeny(int(elem2),Ztmp,Ntmp) #dum2[idx] return string of decay mode, Ztmp, Ntmp, return Z & N of daughter to be ID'd
                            ProgNameL[idx2] = NuclideIdentifier(Ztmp,Ntmp) # ID's daughter name via Ztmp & Ntmp
                    else:
                        DecNameS,Ztmp,Ntmp = translators.Progeny(int(elem1),Ztmp,Ntmp)
                        ProgNameS = NuclideIdentifier(Ztmp,Ntmp)
                DecNameT = str(DecNameS)+str(DecNameL)
                ProgNameT = str(ProgNameS)+str(ProgNameL)
            else:
                DecNameT,Ztmp,Ntmp = translators.Progeny(int(elem),Ztmp,Ntmp)
                ProgNameT = NuclideIdentifier(Ztmp,Ntmp)
            TranslatedMode.append(DecNameT)
            Daughters.append(ProgNameT)
    return(NumOfModes, TranslatedMode, Daughters)

## NEUTRON REACTION SPECIFIC FUNCTIONS
def NeutronRxn_MATID(MatID_file): # obtains MAT ID values for neutron reaction sublibrary
    ID = []
    file1 = open(MatID_file,'r')
    for line in file1:
        if ')' in line:
            a = line.split()
            if len(a[-1]) > 4:
                ID.append(a[-1][-4:])
            else:
                ID.append(a[-1])
    return(ID)

def Get_nRxn(nRxn_file,Z,N):
    nRxnType = []
    nRxnProg = []
    Rxns_not_Tracked = []
    file1 = open(nRxn_file,'r')
    tmp = linecache.getline(nRxn_file,5,None) #pulls line 2 from current_file into cache
    a = tmp.split()
    skip = a[4] # grabs the number of records of descriptive text (NWD in ENDF7.1 manual)
    for i in range(int(skip)+5): #skip over the descriptive records...
        file1.__next__()
    for line in file1: # ...and go straight to the (MF, MT, NC, MOD) descriptive lines
        a = line.split()
        if (a[0]=='3'): # ID's which neutron reactions are tabulated.
            Ztmp,Ntmp,tmp0,tmp1 = translators.MT(int(a[1]),Z,N) # using base Z & N values of isotope, get new Z & N for isotope post neutron reaction
            nRxnType.append(tmp0) # append reaction type
            Rxns_not_Tracked.append(tmp1) # append which reactions are tabulated but not tracked
            nRxnProg.append(NuclideIdentifier(Ztmp,Ntmp)) # append which isotopes correspond to the new Z & N values
        elif ('1  099999' in line): # if you hit the end of the (MF, MT, NC, MOD) records, quit out of function and move onto next nuclide file
            linecache.clearcache()
            return(nRxnType,nRxnProg,Rxns_not_Tracked)

################################################################################
##              START MAIN                                                    ##
################################################################################

#### Check if Isotope List for decay daughter product ID needs to be built.
if not os.path.isfile('IsotopeID.txt'): # if the isotopeID list does not exist, create it
    IsotopeID = open('IsotopeID.txt','w')
    IsotopeID.write('ID \t\t Z, N, \n')
    print('\nBuilding decay sublibrary isotope List.')#Once finished, please re-exeute script and decay library will be built.\n')
    for endf in dec_List:
        dZAID,dID,Z,N = Get_Info(dec_path+'/'+endf)
        IsotopeID.write(str(dID)+'\t\t'+str(Z)+' '+str(N)+' \n')

#### Check if MatID List for neutron RXN needs to be built.
if not os.path.isfile('nRxnMATlist.txt'):
    nRxnMATID = open('nRxnMATlist.txt','w')
    nRxnMATID.write('MAT ID\n')
    print('\nBuidling neutron reaction sublibrary MAT ID list.\n')
    MatID_file = '../ENDF7.1/SubLib_NeutronRxn/n-ENDF-B-VII.1.endf.list'
    MatID = NeutronRxn_MATID(MatID_file)
    for elem in MatID:
        nRxnMATID.write(str(elem)+'\n')

SFyields = os.listdir('../ENDF7.1/SubLib_SFyields/sfy') #gets list of files with Spontaneous fission yields

# Set up xml output information
XML_out = 'DepletionData.xml'
file1 = open(XML_out, 'w', encoding = 'utf-8')
root = ET.Element("Decay_Library", Generator = "INL", Name = "General ENDF7.1", Ver = "1.0")
Lib = ET.ElementTree(root)

#### START LOOPING THROUGH ENDF FILES
dec_count = 0 # start counter for number of files program runs through


dCnt = 0; nCnt=0
while dCnt < len(dec_List):
    decay_filename = dec_List[dCnt]
    print(decay_filename)
    dec_count +=1

    # if (decay_filename == 'dec-002_He_003.endf'):
    #     break

    # decay_file = 'dec-005_B_012.endf'; count = 0

    if ((decay_filename == 'dec-003_Li_008.endf') or (decay_filename == 'dec-003_Li_009.endf')):
        pass
    else:
        decay_file = dec_path+'/'+decay_filename

    ## RADIOACTIVE DECAY INFORMATION.
    dZAID,dID,Z,N = Get_Info(decay_file)
    Half_Life, Mode, Q, BR, nu = Get_DecayInfo(dec_count,decay_file)
    # print(Half_Life, Mode, Q, BR, nu)
    NumOfModes, TranslatedMode, Daughters = TranslateDecayMode(Mode,Z,N,SFyields,decay_file)
    # Daughters = ID_DaughterProducts(Mode,Z,N,SFyields,decay_file)

    ######    PRINT ISOTOPE INFO   #########
    ## RADIOACTIVE DECAY INFORMATION.
    isotope = ET.SubElement(root, ("Isotope"), Halflife = str(Half_Life), Name = str(dID), ZAID = str(dZAID) )
    SubLib_Dec = ET.SubElement(isotope, "Decay")
    Type = ET.SubElement(SubLib_Dec, "type")
    Type.text = str(TranslatedMode).strip("[]")
    Progeny = ET.SubElement(SubLib_Dec, "daughters")
    Progeny.text = str(Daughters).strip("[]")
    Branch_Ratio = ET.SubElement(SubLib_Dec, "branch_ratio")
    Branch_Ratio.text = str(BR).strip("[]")
    ## END RAD DECAY INFO


    while nCnt < len(nRxn_List):
        nRxn_filename = nRxn_List[nCnt]

        ## NEUTRON REACTION INFORMATION.
        nRxn_file = nRxn_path+'/'+nRxn_filename
        nrZAID,nrID,Z,N = Get_Info(nRxn_file)

        if (str(nrID) == str(dID)): #if you have info for both neutron reactions and decay (comparing isotope ID names between decay and neutron reaction sublibraries) do neutron reaction functions
            print('  Match!  '+str(nRxn_filename)+'   adding Neutron Reactions')

            nRxnType, nRxnProg, Rxns_not_Tracked = Get_nRxn(nRxn_file,Z,N)

            ######    PRINT ISOTOPE INFO   #########
            ## NEUTRON REACTION INFORMATION.
            SubLib_nRxn = ET.SubElement(isotope, "Neutron_Reaction")
            nRxnMode = ET.SubElement(SubLib_nRxn, "type")
            nRxnMode.text = str(nRxnType).strip("[]")
            nRxnDau = ET.SubElement(SubLib_nRxn, "daughters")
            nRxnDau.text = str(nRxnProg).strip("[]")
            nRxnNT = ET.SubElement(SubLib_nRxn, "untracked_reactions")
            nRxnNT.text = str(Rxns_not_Tracked).strip("[]")
            ## END NEUTRON REACTION INFO

            nCnt += 1
        else: # if no match between sublibraries, don't do any of the neutron reaction functions and skip to next file in decay sublibrary
            break

    dCnt += 1

Lib.write(file1, encoding='unicode')
file1.close

### temp workaround for parsing.
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
