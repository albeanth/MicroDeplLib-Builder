def Progeny(flag, Z, N):
    operation = {
        0: (0,0),   # gamma
        1: (1,-1),  # beta-
        2: (-1,1),  # e.c./beta+
        3: (0,0),   # I.T.
        4: (-2,-2), # alpha
        5: (0,-1),  # neutron emission
        6: '   SF functionality not available', # Spontaneous fission
        7: (-1,0),  #proton emission
        10: (0,0)   # unknown
    }
    if flag == 6:
        print(operation[flag])
    else:
        Z += operation[flag][0]
        N += operation[flag][1]

    return (Z, N)


def MT(flag, Z, N):
    tmp1 = None
    operation = {
        16: (0,-1,'(n,2n)'),
        17: (0,-2,'(n,3n)'),
        18: ('   fission - functionality not available', '(n,f)'),
        37: (0,-3,'(n,4n)'),
        102: (0, 1,'(n,gamma)'),
        103: (1, -1,'(n,p)'),
        104: (-1, 0,'(n,d)'),
        105: (-1, -1,'(n,t)'),
        106: (-2, 0,'(n,He-3)'),
        107: (-2, -1,'(n,alpha)'),
        108:(-4, -3,'(n,2alpha)'),
        109:(-6, -5,'(n,3alpha)')
    }
    MTlist = operation.keys()
    if flag not in MTlist:
        tmp0 = 'Reaction ' +str(flag)+' not_tracked'; tmp1 = flag
    elif flag == 18:
        print(operation[flag][0])
        tmp0 = operation[flag][1]
    else:
        Z += operation[flag][0]
        N += operation[flag][1]
        tmp0 = operation[flag][2]

    return (Z, N, tmp0, tmp1)
