import sys,os
import linecache
import re

try:
    from colorama import Fore, Back, Style, init
    init(autoreset=True)
    tmpStrY=Fore.YELLOW; tmpStrR=Fore.RED
except ImportError:
    # print('\nYou should get colorama. It\'s pretty sweet.\n')
    tmpStrY=''


def DecProgeny(flag, Z, N,decay_filename): ## Decay mode interpreter - direct from ENDF7.1 manual. Entries 8 and 9 are left out in manual
    FissP_Yield = None
    operation = {
        0: (0, 0, 'gamma'),
        1: (1, -1, 'beta_m'),
        2: (-1, 1, 'e.c./beta_p'),
        3: (0, 0, 'I.T.'),
        4: (-2, -2, 'alpha'),
        5: (0, -1, 'neutron_emission'),
        6: ('spontaneous_fission'),
        7: (-1, 0, 'proton_emission'),
        10: (0, 0, 'unknown'),
    }
    if flag == 6:
        whichLib = 454
        DecName = operation[flag]
        sf_filename = decay_filename.replace('dec-','sfy-')
        sf_path = './ENDF7.1/sfy/'
        sf_file = sf_path+'/'+sf_filename
        sf_List = os.listdir(sf_path)
        if sf_filename not in sf_List:
            print(tmpStrY+'No SF data available for '+str(sf_filename))
            ProgName = 'Unknown, no ENDF distribution given.'
            FissP_Yield = 'Unknown, no ENDF distribution given.'
        else:
            FissP_ID = {}; FissP_Yield = {}
            with open(sf_file,'r') as file1:
                tmp = linecache.getline(sf_file,5,None)
                a = tmp.split()
                skip0 = a[4] #gets number of lines to skip to gloss over "descriptive text"
                tmp = linecache.getline(sf_file,(5+int(skip0)+2),None)
                a = tmp.split()
                skip1 = a[2] # gets number of lines for MT=454 or 459

                for i in range(5+int(skip0)+5): # puts user on line1 of MT=454
                    file1.__next__()

                if whichLib == 459:
                    for i in range(int(skip1)+1): # puts user on line1 of MT=459
                        file1.__next__()
                    ProgName, FissP_Yield = Get_FissionProg(file1,FissP_ID,FissP_Yield)

                elif whichLib == 454:
                    ProgName, FissP_Yield = Get_FissionProg(file1,FissP_ID,FissP_Yield)

                else:
                    print('\n I don\'t know what library you want, please check "whichLib". Quitting.\n')
                    sys.exit()

                ProgName = ProgName['0.000000+0'] #SF data has a "incident" energy of 0. No need to keep key from dictionary. Just pull value.
                FissP_Yield = FissP_Yield['0.000000+0']

    else:
        Z += operation[flag][0]
        N += operation[flag][1]
        DecName = operation[flag][2].strip('')
        ProgName = NuclideIdentifier(Z,N)
    return (DecName, ProgName, FissP_Yield, Z, N)

def MT(flag, Z, N): # ID's neutron reaction types (except fission) and gets daughter products

    ## TO ADD REACTIONS TO BE TRACKED, ADD REACTION TYPE WITH FOLLOWING FORMAT (Z ADJUST ,N ADJUST, ID OF REACTION)
    ## TO REMOVE REACTIONS TO BE TRACKED, JUST COMMENT OUT OR DELETE DICTIONARY ENTRY IN "operation = {}"

    noTrack = ' ' # initialize storing of types that are not being tracked.
    operation = {
        16: (0,-1,'n_2n'),
        17: (0,-2,'n_3n'),
        37: (0,-3,'n_4n'),
        102: (0, 1,'n_gamma'),
        103: (-1, 1,'n_p'),
        104: (-1, 0,'n_d'),
        105: (-1, -1,'n_t'),
        106: (-2, 0,'n_He-3'),
        107: (-2, -1,'n_alpha'),
        108:(-4, -3,'n_2alpha'),
        109:(-6, -5,'n_3alpha')
    }
    MTlist = operation.keys()
    if flag not in MTlist:
        RxnType = ' '; noTrack = flag
    elif flag == 18:
        RxnType = operation[flag][1]
    else:
        Z += operation[flag][0]
        N += operation[flag][1]
        RxnType = operation[flag][2].strip('')

    return (Z, N, RxnType, noTrack)

def MT_fission(flag,nFission_filename): #this ID's fission types, progeny, and yield information for a fissionable isotope
    progeny = []
    CumulYield = []
    ProgList=[]
    operation = {
        18: ('n_f'),
        19: ('n_nf'),
        20: ('n_2nf'),
        21: ('n_3nf'),
        38: ('n_4nf'),
    }
    nFission_Path = './ENDF7.1/nfy'
    nFissionFile = nFission_Path+'/'+nFission_filename
    nFission_List = os.listdir(nFission_Path)
    if nFission_filename not in nFission_List: # not all fissionable isotopes have yield information.
        print(tmpStrY+'No (n,f) data available for '+str(nFission_filename))
        fissType = flag
        FissP_ID = 'Unknown, no ENDF distribution given.'
        FissP_Yield = 'Unknown, no ENDF distribution given.'
    else:
        fissType = operation[flag]
        whichLib = 454 # make either 454 or 459 for independent yield or cumulative yield, respectively.
        with open(nFissionFile,'r') as file1:
            FissP_ID = {}
            FissP_Yield = {}
            tmp = linecache.getline(nFissionFile,5,None)
            a = tmp.split()
            skip0 = a[4] #gets number of lines to skip to gloss over "descriptive text"
            tmp = linecache.getline(nFissionFile,(5+int(skip0)+2),None)
            a = tmp.split()
            skip1 = a[2] # gets number of lines for MT=454 or 459

            for i in range(5+int(skip0)+5): # puts user on line1 of MT=454
                file1.__next__()

            if whichLib == 459:
                for i in range(int(skip1)+1): # puts user on line1 of MT=459
                    file1.__next__()
                FissP_ID,FissP_Yield = Get_FissionProg(file1,FissP_ID,FissP_Yield)

            elif whichLib == 454:
                FissP_ID,FissP_Yield = Get_FissionProg(file1,FissP_ID,FissP_Yield)

            else:
                print(tmpStrG+'\n I don\'t know what library you want, please check "whichLib". Quitting.\n')
                sys.exit()

    return (fissType, FissP_ID, FissP_Yield)

def Get_FissionProg(file1,FissP_ID,FissP_Yield): # takes either spontaneous fission or neutron induced fission and ID's fission products and their yields
    line = file1.readline() # read line 1 of MT=454/459
    a = line.split()
    NumEnergySets = a[2] # ID how many energy dependent fission product yields are given -- "LE" in ENDF7.1 manual, Ch.8
    line = file1.readline() # read line 2 of MT=454/459
    EScount = 0 # energy set counter
    while EScount<int(NumEnergySets):
        ProgList = [] # re-initialize progeny list for each energy set
        a = line.split()
        Energy = a[0] # grab incident energy...
        length = a[4] # ... and number of expected entries.
        for line in file1: # start looping over incident energy yield data set.
            if len(ProgList) >= int(length):
                break
            tmp = line.split() # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+09025',  '8454', '3']
            tmp = tmp[:-2]     # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+09025']
            if len(tmp) < 6: # can be at the end of energy set if list formatting ends before line ends
                del tmp[-1]
            else:
                tmp.append(tmp[-1][:-4]) # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+09025', '0.000000+0']
                del tmp[5]        # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+0']
            ProgList.extend(tmp)

        ID = []; Yield = []; StDes = []
        i=0; j=1; k=2 # i = index for nuclide identifier; j = index for state designator; k = index for cumulative yield of isotope i
        while(i<len(ProgList)):
            pdDm = ProgList[i] # ZAFP of fission product
            sdDum = ProgList[j] #state designator
            tmp2 = NuclideIdentifier_fission(pdDm)
            if sdDum != '0.000000+0': #metastable state.
                tmp2 = str(tmp2)+'m'
            else:
                pass
            ID.append(tmp2)
            cyDum = ProgList[k]
            Yield.append(ScientificNotation(cyDum))
            i+=4; j+=4; k+=4 # see C_n(E_i) in Ch8 of ENDF7.1 manual for description

        FissP_ID[Energy] = ID
        FissP_Yield[Energy] = Yield

        EScount+=1

    return(FissP_ID,FissP_Yield)

def NuclideIdentifier_fission(ZAFP): # takes ZAFP value from fission product ID (MT=454/459) and identifies the daughter product explicitely.
    with open('IsotopeID.txt','r') as IsotopeList:
        for line in IsotopeList:
            a = line.split()
            if (str(ZAFP)==a[3]):
                name=a[0]
                # print('Fission Product = '+str(name))
                break
            else:
                name = str(ZAFP)
                ## uncomment code below and comment single line about to print unknown fission products in standard ZAID form
                # ZA = re.split("([\+\-])", ZAFP) #splits ZA into 3 components
                # ZA.insert(1, "E") #inserts the E needed for scientific notation
                # ZA = "".join(ZA[0:]) #rejoins array of strings into single string
                # name = int(float(ZA)*10.0)
    return(name)# once the daughter is found, break out or else it will cycle through the entire file every single time is searches

def NuclideIdentifier(Z,N): # takes adjust Z & N values from decay type and identifies the daughter product explicitely.
    with open('IsotopeID.txt','r') as IsotopeList:
        for line in IsotopeList:
            a = line.split()
            if ((str(Z)==a[1]) and (str(N)==a[2])):
                name=a[0]
                break # once the daughter is found, break out or else it will cycle through the entire file every single time is searches and waste time
            elif (((Z==3) and (N==5)) or ((Z==3) and (N==6))): #if Li-8 or Li-9 (wrong file format, ticket submitted)
                name = 'N/A'
    return(name)

def ScientificNotation(tmp): # convert ENDF "1.00000+2" to "1.00000E+2"
    if (('+' in tmp) or ('-' in tmp)):
        tmp = re.split("([\+\-])", tmp) #splits tmp into 3 components
        tmp.insert(1, "E") #inserts the "E" needed for scientific notation
        tmp = "".join(tmp[0:]) #rejoins array of strings into single string
    else:
        pass
    return(tmp)
