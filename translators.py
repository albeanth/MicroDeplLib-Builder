import sys,os
import linecache
import re

def Progeny(flag, Z, N): ## Decay mode interpreter - direct from ENDF7.1 manual. Entries 8 and 9 are left out in manual
    operation = {
        0: (0, 0, 'gamma'),
        1: (1, -1, 'beta-'),
        2: (-1, 1, 'e.c./beta+'),
        3: (0, 0, 'I.T.'),
        4: (-2, -2, 'alpha'),
        5: (0, -1, 'neutron emission'),
        6: ('   SF functionality not available', 'Spontaneous fission'),
        7: (-1,0, 'proton emission'),
        10: (0,0, 'unknown'),
    }
    if flag == 6:
        print(operation[flag][0])
        name = operation[flag][1]
    else:
        Z += operation[flag][0]
        N += operation[flag][1]
        name = operation[flag][2]

    return (name, Z, N)


def MT(flag, Z, N):
    tmp1 = ' '
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
        # tmp0 = 'Reaction ' +str(flag)+' not tracked'; tmp1 = flag
        tmp0 = ' '; tmp1 = flag
    elif flag == 18:

        print(operation[flag][0])
        tmp0 = operation[flag][1]
    else:
        Z += operation[flag][0]
        N += operation[flag][1]
        tmp0 = operation[flag][2]

    return (Z, N, tmp0, tmp1)

def MT_fission(flag,nFission_filename):
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
    nFission_Path = '../ENDF7.1/SubLib_nfYields/nfy'
    nFissionFile = nFission_Path+'/'+nFission_filename
    nFission_List = os.listdir(nFission_Path)
    if nFission_filename not in nFission_List:
        fissType = flag
        FissP_ID = {'Unknown': 'Unknown. No ENDF distribution given.'}
        FissP_Yield = {'Unknown': 'Unknown. No ENDF distribution given.'}
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

            for i in range(5+int(skip0)+5): # puts user on line2 of MT=454
                file1.__next__()

            if whichLib == 459:
                for i in range(int(skip1)+1): # puts user on line1 of MT=459
                    file1.__next__()
                FissP_ID,FissP_Yield = Get_FissionProg(file1,FissP_ID,FissP_Yield)

            elif whichLib == 454:
                FissP_ID,FissP_Yield = Get_FissionProg(file1,FissP_ID,FissP_Yield)

    return (fissType, FissP_ID, FissP_Yield)

def Get_FissionProg(file1,FissP_ID,FissP_Yield):
    line = file1.readline() # read line 1 and grab number of yield sets
    a = line.split()
    NumEnergySets = a[2]
    line = file1.readline() # read line 2...

    EScount = 0
    while EScount<int(NumEnergySets):
        ProgList = []
        a = line.split()
        Energy = a[0] # ... and grab incident energy...
        length = a[4] # ... and number of expected entries.
        for line in file1: # start looping over incident energy yield data. should produce ProgList of length == "length = a[4]"
            if len(ProgList) >= int(length):
                break
            tmp = line.split() # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+09025',  '8454', '3']
            tmp1=tmp
            tmp = tmp[:-2]     # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+09025']
            if len(tmp) < 6:
                del tmp[-1]
            else:
                tmp.append(tmp[-1][:-4]) # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+09025', '0.000000+0']
                del tmp[5]        # ['2.406600+4', '0.000000+0', '0.000000+0', '0.000000+0', '2.406700+4', '0.000000+0']
            ProgList.extend(tmp)

        # Progeny[Energy] = ProgList

        linecache.clearcache()

        ID = []
        Yield = []
        i=0; j=1; k=2 # i = index for nuclide identifier; j = index for state designator; k = index for cumulative yield of isotope i
        while(i<len(ProgList)):
            pdDm = ProgList[i]
            # print(pdDm)
            ID.append(NuclideIdentifier_fission(pdDm))
            sdDum = ProgList[j]
            cyDum = ProgList[k]
            Yield.append(ScientificNotation(cyDum))
            i+=4; j+=4; k+=4

        FissP_ID[Energy] = ID
        FissP_Yield[Energy] = Yield

        EScount+=1

    return(FissP_ID,FissP_Yield)

def NuclideIdentifier_fission(ZAFP):
    # print(ZAFP)
    with open('IsotopeID.txt','r') as IsotopeList:
        for line in IsotopeList:
            a = line.split()
            if (str(ZAFP)==a[3]):
                name=a[0]
                # print('Fission Product = '+str(name))
                break
            else:
                name = str(ZAFP)
                # ZA = re.split("([\+\-])", ZAFP) #splits ZA into 3 components
                # ZA.insert(1, "E") #inserts the E needed for scientific notation
                # ZA = "".join(ZA[0:]) #rejoins array of strings into single string
                # name = int(float(ZA)*10.0)
    return(name)# once the daughter is found, break out or else it will cycle through the entire file every single time is searches



def NuclideIdentifier(Z,N):
    with open('IsotopeID.txt','r') as IsotopeList:
        for line in IsotopeList:
            a = line.split()
            if ((str(Z)==a[1]) and (str(N)==a[2])):
                name=a[0]
                break # once the daughter is found, break out or else it will cycle through the entire file every single time is searches
            elif (((Z==3) and (N==5)) or ((Z==3) and (N==6))): #if Li-8 or Li-9 (wrong file format, ticket submitted)
                name = 'N/A'
    return(name)

def ScientificNotation(tmp): # convert ENDF "1.00000+2" to "1.00000E+2"
    if (('+' in tmp) or ('-' in tmp)):
        tmp = re.split("([\+\-])", tmp) #splits ModeTmp into 3 components
        tmp.insert(1, "E") #inserts the E needed for scientific notation
        tmp = "".join(tmp[0:]) #rejoins array of strings into single string
    else:
        pass
    return(tmp)
