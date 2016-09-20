import sys,os
import os.path
import linecache
import re
import numpy as np

try:
    from colorama import Fore, Back, Style, init
    init(autoreset=True)
    tmpStrY=Fore.YELLOW; tmpStrR=Fore.RED; tmpStrG=Fore.GREEN
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    tmpStrY=''

#packages needed for xml creation and parsing
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

# accompanying python file for required functions
import translators as trls

import zipfile
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

## GENERAL PURPOSE FUCNTIONS
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
    AWR = trls.ScientificNotation(AWRtmp)
    A = float(AWR) * mass_n # use to get Atomic Mass
    Z = round((float(ZA) - A)/1000.) # calculates the atomic number
    N = round(A-Z) # calculates the number of neutrons
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

    return (ZAtmp,ZAID,ID,Z,N)

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
                print(tmpStrR+'it''s a polynomial? really? ...quitting.')
                sys.exit()
            elif LNU ==2:
                tmp = linecache.getline(current_file,(nuLibStart+2))
                a = tmp.split()
                nu = a[0]

        elif " 8457" in line: #signifies start of library 8 (decay) MT=457, rad decay data set
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
                Half_Life=trls.ScientificNotation(a[0])
                # get number of decay modes (NDK)
                tmp = linecache.getline(current_file,(DecLibStart+3),None)
                a = tmp.split()
                NDK=a[5]
                if ((count>999) & (NDK.endswith(str(count)))):
                    NDK = int(NDK[:-len(str(count))])
                else:
                    NDK=int(NDK)
                # get decay info: Modes, Q values, and Branching Ratios (BR)
                Mode = np.empty(NDK,dtype=float) #create decay mode vector of length equal to NDK
                Q = np.empty(NDK,dtype=float) #create Q value vector for each mode
                BR = np.empty(NDK,dtype=float) #create Branching ratio  vector of length equal to NDK
                RTYP_start = DecLibStart+4 # line number in which mode of decay identifiers (RYTP) start
                i = 0
                while (i<NDK):
                    tmp = linecache.getline(current_file,(RTYP_start+i),None)
                    a = tmp.split()

                    ModeTmp = a[0]
                    Mode[i] = trls.ScientificNotation(ModeTmp)

                    QTmp = a[2]
                    Q[i] = trls.ScientificNotation(QTmp)

                    BRTmp = a[4]
                    BR[i] = trls.ScientificNotation(BRTmp)

                    i += 1
                break
    return(Half_Life, Mode, Q, BR, nu)

def TranslateDecayMode(Mode,Z,N,current_file,decay_filename): #translate RTYP numbers to readable decay types
    if Mode == None:
        NumOfModes = None
        Daughters = None
        sfYield = None
    else:
        NumOfModes=len(Mode)
        Daughters = {}
        sfYield = []
        for idy,elem in enumerate(Mode):
            Ztmp = Z; Ntmp = N
            # not all decay types for all isotopes lead to actual results. the following if...elif... statements ID which ones are issues.
            if ((current_file == './ENDF7.1/decay/dec-028_Ni_048.endf') and (elem == 2.0)):
                print(tmpStrY+'  No data for Ni-48 beta+ decay. Need data for Co-48. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-028_Ni_048.endf') and (elem == 7.7)):
                print(tmpStrY+'  No data for proton-proton decay. Need data for Co-47. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-098_Cf_239.endf') and (elem == 2.0)):
                print(tmpStrY+'  No data for Bk-239 from Cf-239 beta+ decay. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-098_Cf_256.endf') and (elem == 4.0)):
                print(tmpStrY+'  No data for Cf-256 alpha decay (and it\'s probability is 1E-9). Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-099_Es_240.endf') and (elem == 4.0)):
                print(tmpStrY+'  No data for Es-240 alpha decay. Need data for Bk-236. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-099_Es_243.endf') and (elem == 4.0)):
                print(tmpStrY+'  No data for Es-243 alpha decay. Need data for Bk-239. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-099_Es_258.endf') and (elem == 2.0)):
                print(tmpStrY+'  No data for Es-243 Beta+ decay. Need data for Cf-258 (which doesn\'t exist?). Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-104_Rf_253.endf') and (elem == 4.0)):
                print(tmpStrY+'  No data for Rf-253 alpha decay. Need data for No-249. Passing...')
                tmp = 'N/A'
                continue
            elif ((current_file == './ENDF7.1/decay/dec-110_Ds_279m1.endf') and (elem == 4.0)):
                print(tmpStrY+'  No data for Ds-279m1 alpha decay. Need data for Hs-275. Passing...')
                tmp = 'N/A'
                continue
            if ((elem-int(elem))>0.): # checks to see if there are multiple decays, e.g. 1.534
                tmp1 = str(elem).split('.') # splits decimal number into two strings, e.g. ['1', '534']
                tmp2 = list(map(str, tmp1[1])) # splits up '534' into independent strings - e.g. ['5', '3', '4']
                tmp1[1]=tmp2 # sets tmp2 into tmp1[1] --> tmp1 = ['1', ['5', '3', '4']]
                for idx,elem1 in enumerate(tmp1):
                    if elem1.count(elem1)==0: # if elem1 is the ['5', '3', '4'] term
                        DecNameL = np.chararray(len(elem1),16,True)
                        ProgNameL = np.chararray(len(elem1),16,True)
                        for idx2,elem2 in enumerate(elem1): # for each term, '5','3','4'
                            DecNameL[idx2],ProgNameL[idx2],dummy, Ztmp, Ntmp = trls.DecProgeny(int(elem2),Ztmp,Ntmp,decay_filename)
                    else:
                        DecNameS,ProgNameS,dummy, Ztmp, Ntmp = trls.DecProgeny(int(elem1),Ztmp,Ntmp,decay_filename)
                DecNameT = str(DecNameS)+str(DecNameL) # 'S' -> short (1.) 'L' -> Long (.534) therefore, 'T' -> total (1.534)
                ProgNameT = str(ProgNameS)+str(ProgNameL)
            else:
                DecNameT,ProgNameT,dummy, Ztmp, Ntmp = trls.DecProgeny(int(elem),Ztmp,Ntmp,decay_filename)
            Daughters[DecNameT] = ProgNameT
            sfYield = dummy
    return(NumOfModes, Daughters, sfYield)

## NEUTRON REACTION SPECIFIC FUNCTIONS
def Get_nRxn(nRxn_file,Z,N,nFission_filename): # this goes through each neutron reaction file, ID's reaction types, and daughter products. nRxnTypes, fission products and yields returned as dictionaries and nRxns not tracked are returned as a list.
    nRxnType = {}
    FissP_ID = {}
    FissP_Yield = {}
    Rxns_not_Tracked = []
    with open(nRxn_file,'r') as file1:
        tmp = linecache.getline(nRxn_file,5,None) #pulls line 2 from current_file into cache
        a = tmp.split()
        skip = a[4] # grabs the number of records of descriptive text (NWD in ENDF7.1 manual)
        for i in range(int(skip)+5): #skip over the descriptive records...
            file1.__next__()
        for line in file1: # ...and go straight to the (MF, MT, NC, MOD) descriptive lines
            a = line.split()
            if (a[0]=='3'): # ID's which neutron reactions are tabulated.
                if (a[1]=='18'): #or (a[1]=='19') or (a[1]=='20') or (a[1]=='21') or (a[1]=='38')
                    fissType, FissP_ID, FissP_Yield = trls.MT_fission(int(a[1]),nFission_filename)
                else:
                    Ztmp,Ntmp,RxnType,noTrack = trls.MT(int(a[1]),Z,N) # using base Z & N values of isotope, get new Z & N for isotope post neutron reaction
                    if noTrack == ' ':
                        pass
                    else:
                        Rxns_not_Tracked.append(noTrack) # append which reactions are tabulated but not tracked
                    if RxnType == ' ':
                        pass
                    else:
                        nRxnType[RxnType] = trls.NuclideIdentifier(Ztmp,Ntmp)

            elif ('1  099999' in line): # if you hit the end of the (MF, MT, NC, MOD) records, quit out of function and move onto next nuclide file
                linecache.clearcache()
                return(nRxnType,FissP_ID, FissP_Yield,Rxns_not_Tracked)


################################################################################
##              START MAIN                                                    ##
################################################################################

## SET PATHS OF SUBLIBRARIES
dec_path = './ENDF7.1/decay' # set path to point to ENDF7.1 files
nRxn_path = './ENDF7.1/neutrons' # set path to point to ENDF7.1 files
# be sure to adjust path for neutron induced fission in MT_fission() in translators.py

## DEFINE CONSTANTS
mass_n = 1.00866491578 # mass of nuetron in amu. Source - ENDF7.1 manual

try: # if this throws an error than make the directory and create the files
    os.listdir('./ENDF7.1')
except FileNotFoundError:
    os.mkdir('ENDF7.1')

if len(os.listdir('./ENDF7.1')) == 0: # if the ENDF7.1 directory is empty
    print('Obataining Needed NNDC ENDF7.1 Sublibraries')
    cwd = os.getcwd() # return string of current working directory
    # print(cwd)
    sys.path.insert(0, os.path.join(cwd, '..'))

    baseUrl = 'http://www.nndc.bnl.gov/endf/b7.1/zips/'

    files = ['ENDF-B-VII.1-neutrons.zip',
             'ENDF-B-VII.1-decay.zip',
             'ENDF-B-VII.1-sfy.zip',
             'ENDF-B-VII.1-nfy.zip']
    block_size = 1024 # specifies the chunks in which the progress bar prints itself. can change around...

    filesComplete = []
    for f in files:
        # Establish connection to URL
        url = baseUrl + f
        req = urlopen(url)

        # Get file size from header
        if sys.version_info[0] < 3: #check to see what version of python you have
            file_size = int(req.info().getheaders('Content-Length')[0])
        else:
            file_size = req.length
            # print('size = '+str(file_size))
        downloaded = 0

        # Check if file already downloaded
        if os.path.exists(f):
            if os.path.getsize(f) == file_size:
                print('Skipping... Already downloaded ' + f)
                filesComplete.append(f)
                continue
            else:
                if sys.version_info[0] < 3:
                    overwrite = raw_input('Overwrite {0}? ([y]/n) '.format(f))
                else:
                    overwrite = input('Overwrite {0}? ([y]/n) '.format(f))
                if overwrite.lower().startswith('n'):
                    continue

        # Copy file to disk
        print('Downloading {0}... '.format(f), end='')
        with open(f, 'wb') as fh:
            while True:
                chunk = req.read(block_size)
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                status = '{0:10}  [{1:3.2f}%]'.format(downloaded, downloaded * 100. / file_size)
                print(status + chr(8)*len(status), end='')
            print('')
            filesComplete.append(f)

        if f not in filesComplete:
            continue

    for f in filesComplete:
        print('...extracting '+str(f))
        zip_ref = zipfile.ZipFile(str(cwd)+'/'+str(f), 'r')
        zip_ref.extractall(path=str(cwd)+'/ENDF7.1')
        zip_ref.close()
        print('removing .zip file -> '+str(f))
        os.remove(f)
else:
    pass

try:
    dec_List = os.listdir(dec_path) # get list of isotope files in decay sublibrary
    nRxn_List = os.listdir(nRxn_path) # get list of isotope files in neutron reactions
except FileNotFoundError:
    pass

#### Check if Isotope List for decay daughter product ID needs to be built.
if not os.path.isfile('IsotopeID.txt'): # if the isotopeID list does not exist, create it
    IsotopeID = open('IsotopeID.txt','w')
    IsotopeID.write('ID \t\tZ   N   ZA\n')
    print(tmpStrG + '\nBuilding isotope List from decay sublibrary.')#Once finished, please re-exeute script and decay library will be built.\n')
    for endf in dec_List:
        ZAtmp, dZAID,dID,Z,N = Get_Info(dec_path+'/'+endf)
        IsotopeID.write(str(dID)+'\t\t'+str(Z)+'   '+str(N)+'   '+str(ZAtmp)+'\n')
    IsotopeID.close()


# Set up xml output information
XML_out = 'DepletionData.xml'
file1 = open(XML_out, 'w', encoding = 'utf-8')
root = ET.Element("Decay_Library", Generator = "INL", Name = "General ENDF7.1", Ver = "1.0")
Lib = ET.ElementTree(root)



#### START LOOPING THROUGH ENDF FILES
ParamCnt = 0
dCnt = 0; nCnt=0
while dCnt < len(dec_List):
    decay_filename = dec_List[dCnt]
    dCnt += 1

    ## uncomment following two lines to stop library building with the listed filename
    # if (decay_filename == 'dec-064_Gd_147.endf'):
    #     break

    print(Style.DIM + decay_filename+'   #'+str(dCnt))

    if ((decay_filename == 'dec-003_Li_008.endf') or (decay_filename == 'dec-003_Li_009.endf')):
        continue
    else:
        decay_file = dec_path+'/'+decay_filename

    ################  RADIOACTIVE DECAY INFORMATION  ################
    ZAtmp,dZAID,dID,Z,N = Get_Info(decay_file)
    Half_Life, Mode, Q, BR, nu = Get_DecayInfo(dCnt,decay_file)
    NumOfModes, Daughters, sfYield = TranslateDecayMode(Mode,Z,N,decay_file,decay_filename)

    ## PRINT RADIOACTIVE DECAY INFORMATION
    isotope = ET.SubElement(root, ("Isotope"), Decay_Constant = math.log(2)/str(Half_Life), Name = str(dID), ZAID = str(dZAID) )
    if Daughters == None:
        pass
    else:
        SubLib_Dec = ET.SubElement(isotope, "Decay")
        for idx,pair in enumerate(sorted(Daughters.items())):
            ParamCnt += 1
            tmpstr = 'decMode_'+str(idx)
            tmpstr = ET.SubElement(SubLib_Dec, "Mode_"+str(idx+1))
            tmpstr.attrib['AType'] = str(pair[0])
            tmpstr.attrib['Branch_Ratio'] = str(BR[idx]).strip("[]")
            tmpstr.attrib['Daughters'] = str(pair[1])
            if not sfYield==None:
                tmpstr.attrib['SF_Yield'] = str(sfYield).strip("[]")
            tmpstr.attrib['Q_value'] = str(Q[idx]).strip("[]")
    ################  END RAD DECAY INFO  ################

    ################  NEUTRON REACTION INFORMATION  ################
    while nCnt < len(nRxn_List):
        nRxn_filename = nRxn_List[nCnt]
        nFission_filename = nRxn_filename.replace('n-','nfy-')

        nRxn_file = nRxn_path+'/'+nRxn_filename
        ZAtmp,nrZAID,nrID,Z,N = Get_Info(nRxn_file)

        if (str(nrID) == str(dID)): #if you have info for both neutron reactions and decay (comparing isotope ID names between decay and neutron reaction sublibraries) do neutron reaction functions
            print(Style.DIM + '  ....adding Neutron Reactions!'+' (file #'+str(nCnt)+')')
            nCnt += 1
            nRxnType, FissProg, FissYield, Rxns_not_Tracked = Get_nRxn(nRxn_file,Z,N,nFission_filename)

            ## NEUTRON REACTION INFORMATION
            SubLib_nRxn = ET.SubElement(isotope, "Neutron_Reaction")
            for idx,pair in enumerate(sorted(nRxnType.items())):
                ParamCnt +=1
                tmpstr = 'nRxnMode_'+str(idx)
                tmpstr = ET.SubElement(SubLib_nRxn, str(pair[0]))
                tmpstr.text = str(pair[1])

            listing = sorted(FissProg.keys())
            if len(listing) > 0:
                nFissDau = ET.SubElement(SubLib_nRxn, "fission_daughters")
                nFissDau.text = str(FissProg[listing[0]])
                nFissY = ET.SubElement(SubLib_nRxn, "fission_yield")
                nFissY.text = str(FissYield[listing[0]])
            else:
                pass

            nRxnNT = ET.SubElement(SubLib_nRxn, "untracked_reactions")
            nRxnNT.text = str(Rxns_not_Tracked).strip("[]")

        else: # if no match between sublibraries, don't do any of the neutron reaction functions and skip to next file in decay sublibrary
            break

    ################  END NEUTRON REACTION INFO  ##################


Lib.write(file1, encoding='unicode')
file1.close

print(tmpStrG+'\nDone building library.\n')

# print('Number of Parameters in Library = '+str(ParamCnt))

### temp workaround for parsing.
# os.system('xmllint --format DecayData.xml > DecayData_tmp.xml')
# os.system('rm DecayData.xml')
# os.system('mv DecayData_tmp.xml DecayData.xml')

# doc2 = minidom.parse(XML_out)

#### The code below works on simple toy problems but not here.... so use workaround system commands above.
## use xml.dom.minidom and parse out previously created xml file.
# xml = minidom.parse('DepletionData.xml')
# pretty_xml_as_string = xml.toprettyxml()
# file2 = open('DepletionData.xml', 'w', encoding = 'utf-8')
# file2.write(pretty_xml_as_string)
# file2.close()
