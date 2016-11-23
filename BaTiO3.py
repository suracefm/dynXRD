
def structure_BaTiO3():
    cryst = {} # cryst.copy()

    cryst['name']   = b"BaTiO3"
    cryst['a']      =  3.9999
    cryst['b']      =  3.9999
    cryst['c']      =  4.0170
    cryst['alpha']  = 90.0
    cryst['beta']   = 90.0
    cryst['gamma']  = 90.0
    cryst['n_atom'] = 5
    cryst['atom']   = [{'x': 0.0,  'y': 0.0,  'z': 0.0,  'Zatom': 56, 'fraction': 1.0, 'label':'Ba'},
                       {'x': 0.5,  'y': 0.5,  'z': 0.4820,  'Zatom':22, 'fraction': 1.0, 'label':'Ti'},
                       {'x': 0.5,  'y': 0.5,  'z': 0.0160,  'Zatom': 8, 'fraction': 1.0, 'label':'O'},
                       {'x': 0.5,  'y': 0.0,  'z': 0.5150,  'Zatom': 8, 'fraction': 1.0, 'label':'O'},
                       {'x': 0.0,  'y': 0.5,  'z': 0.5150,  'Zatom': 8, 'fraction': 1.0, 'label':'O'}]
    cryst['volume'] = 64.269

    return cryst
