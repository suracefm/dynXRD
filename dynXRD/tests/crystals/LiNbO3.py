def structure_LiNbO3():
    cryst = {} # cryst.copy()

    cryst['name']   = b"LiNbO3"
    cryst['a']      =  5.14829
    cryst['b']      =  5.14829
    cryst['c']      =  13.8631
    cryst['alpha']  = 90.0
    cryst['beta']   = 90.0
    cryst['gamma']  = 120.0
    cryst['initialize_atom']   = [{'x': 0.0,  'y': 0.0,  'z': 0.2829,  'Zatom': 3, 'fraction': 1.0, 'label':'Li'},
                       {'x': 0.0,  'y': 0.0,  'z': 0.0,  'Zatom':41, 'fraction': 1.0, 'label':'Nb'},
                       {'x': 0.0492,  'y': 0.3446,  'z': 0.0647,  'Zatom': 8, 'fraction': 1.0, 'label':'O'},
                       ]

    symmetries=[['-x+y', 'y', 'z+1/2'],['x', 'x-y', 'z+1/2'],['-y', '-x', 'z+1/2'],['-x+y', '-x', 'z'],['-y', 'x-y', 'z'],\
               ['x', 'y', 'z'],['-x+y+2/3', 'y+1/3', 'z+5/6'],['x+2/3', 'x-y+1/3', 'z+5/6'],['-y+2/3', '-x+1/3', 'z+5/6'],\
               ['-x+y+2/3','-x+1/3', 'z+1/3'], ['-y+2/3', 'x-y+1/3', 'z+1/3'],['x+2/3', 'y+1/3', 'z+1/3'],\
               ['-x+y+1/3', 'y+2/3', 'z+1/6'],['x+1/3', 'x-y+2/3', 'z+1/6'], ['-y+1/3', '-x+2/3', 'z+1/6'],\
               ['-x+y+1/3', '-x+2/3', 'z+2/3'],['-y+1/3', 'x-y+2/3', 'z+2/3'],['x+1/3', 'y+2/3', 'z+2/3']]

    atomset=set([])
    for dict in cryst['initialize_atom']:
        x, y, z, Zatom, fraction, label=(dict.get(key) for key in ['x', 'y', 'z', 'Zatom', 'fraction', 'label'])
        for expr in symmetries:
            atomset.add((eval(expr[0]),eval(expr[1]),eval(expr[2]),Zatom, fraction, label))

    cryst['atom']=[{'x':atom[0], 'y':atom[1], 'z':atom[2], 'Zatom':atom[3], 'fraction':atom[4], 'label':atom[5]} for atom in atomset]
    cryst['n_atom'] =  len(cryst['atom'])


    cryst['volume'] = 318.21

    return cryst
