import numpy as np

def make_h_param_vector(twoq_inf, oneq_inf, idle_inf):
    parameters = {}
    
    parameter_indexing = []
    for gatename in ['Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2']:
        parameters[gatename] = {}
        parameters[gatename][('H', 'X')] = oneq_inf/3
        parameters[gatename][('H', 'Z')] = oneq_inf/3
        parameters[gatename][('H', 'Y')] = oneq_inf/3
        parameter_indexing += [(gatename,'X'),(gatename,'Y'),(gatename,'Z')]

    parameters['Gi'] = {}
    parameters['Gi'][('H', 'X')] = oneq_inf/3
    parameters['Gi'][('H', 'Z')] = oneq_inf/3
    parameters['Gi'][('H', 'Y')] = oneq_inf/3
    parameter_indexing += [('Gi','X'),('Gi','Y'),('Gi','Z')]    

    parameters['Gc0'] = {}
    parameters['Gc0'][('H', 'X')] = idle_inf/3
    parameters['Gc0'][('H', 'Z')] = idle_inf/3
    parameters['Gc0'][('H', 'Y')] = idle_inf/3
    parameter_indexing += [('Gc0','X'),('Gc0','Y'),('Gc0','Z')]   
    
    parameters['Gcphase'] = {}
    parameters['Gcphase'][('H', 'IX')] = twoq_inf/15
    parameters['Gcphase'][('H', 'IY')] = twoq_inf/15
    parameters['Gcphase'][('H', 'IZ')] = twoq_inf/15
    
    parameters['Gcphase'][('H', 'XI')] = twoq_inf/15
    parameters['Gcphase'][('H', 'XX')] = twoq_inf/15
    parameters['Gcphase'][('H', 'XY')] = twoq_inf/15
    parameters['Gcphase'][('H', 'XZ')] = twoq_inf/15
    
    parameters['Gcphase'][('H', 'YI')] = twoq_inf/15
    parameters['Gcphase'][('H', 'YX')] = twoq_inf/15
    parameters['Gcphase'][('H', 'YY')] = twoq_inf/15
    parameters['Gcphase'][('H', 'YZ')] = twoq_inf/15
    
    parameters['Gcphase'][('H', 'ZI')] = twoq_inf/15
    parameters['Gcphase'][('H', 'ZX')] = twoq_inf/15
    parameters['Gcphase'][('H', 'ZY')] = twoq_inf/15
    parameters['Gcphase'][('H', 'ZZ')] = twoq_inf/15

    #parameters['Mdefault'] = {('S', 'X'): spam_error}
    #parameters['rho0'] = {('S', 'X'): spam_error}
    
    
    # Parameter indexing, used throughout the analysis. It defines our "theta" vector
   
    return parameters
    
def make_depol_error_param_vector(twoq_inf, oneq_inf, spam_error, idle_inf):
    parameters = {}
    
    parameter_indexing = []
    for gatename in ['Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2']:
        parameters[gatename] = {}
        parameters[gatename][('S', 'X')] = oneq_inf/3
        parameters[gatename][('S', 'Z')] = oneq_inf/3
        parameters[gatename][('S', 'Y')] = oneq_inf/3
        parameter_indexing += [(gatename,'X'),(gatename,'Y'),(gatename,'Z')]

    parameters['Gi'] = {}
    parameters['Gi'][('S', 'X')] = oneq_inf/3
    parameters['Gi'][('S', 'Z')] = oneq_inf/3
    parameters['Gi'][('S', 'Y')] = oneq_inf/3
    parameter_indexing += [('Gi','X'),('Gi','Y'),('Gi','Z')]    

    parameters['Gc0'] = {}
    parameters['Gc0'][('S', 'X')] = oneq_inf/3
    parameters['Gc0'][('S', 'Z')] = oneq_inf/3
    parameters['Gc0'][('S', 'Y')] = oneq_inf/3
    parameter_indexing += [('Gc0','X'),('Gc0','Y'),('Gc0','Z')]   
    
    parameters['Gcphase'] = {}
    parameters['Gcphase'][('S', 'IX')] = twoq_inf/15
    parameters['Gcphase'][('S', 'IY')] = twoq_inf/15
    parameters['Gcphase'][('S', 'IZ')] = twoq_inf/15
    
    parameters['Gcphase'][('S', 'XI')] = twoq_inf/15
    parameters['Gcphase'][('S', 'XX')] = twoq_inf/15
    parameters['Gcphase'][('S', 'XY')] = twoq_inf/15
    parameters['Gcphase'][('S', 'XZ')] = twoq_inf/15
    
    parameters['Gcphase'][('S', 'YI')] = twoq_inf/15
    parameters['Gcphase'][('S', 'YX')] = twoq_inf/15
    parameters['Gcphase'][('S', 'YY')] = twoq_inf/15
    parameters['Gcphase'][('S', 'YZ')] = twoq_inf/15
    
    parameters['Gcphase'][('S', 'ZI')] = twoq_inf/15
    parameters['Gcphase'][('S', 'ZX')] = twoq_inf/15
    parameters['Gcphase'][('S', 'ZY')] = twoq_inf/15
    parameters['Gcphase'][('S', 'ZZ')] = twoq_inf/15

    parameters['Mdefault'] = {('S', 'X'): spam_error}
    parameters['rho0'] = {('S', 'X'): spam_error}
    
    
    # Parameter indexing, used throughout the analysis. It defines our "theta" vector
   
    return parameters
    
def make_twirled_error_param_vector(scale=1):
    parameters = {}
    
    parameter_indexing = []
    for gatename in ['Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2',]:
        parameters[gatename] = {}
        parameters[gatename][('S', 'X')] = 0.2e-3*scale
        parameters[gatename][('S', 'Z')] = 0.6e-3*scale 
        parameters[gatename][('S', 'Y')] = 0.2e-3*scale
        parameter_indexing += [(gatename,'X'),(gatename,'Y'),(gatename,'Z')]
    
    parameters['Gcphase'] = {}
    parameters['Gcphase'][('S', 'IX')] = 5e-5*scale
    parameters['Gcphase'][('S', 'IY')] = 5e-5*scale
    parameters['Gcphase'][('S', 'IZ')] = 5e-4*scale
    
    parameters['Gcphase'][('S', 'XI')] = 5e-5*scale
    parameters['Gcphase'][('S', 'XX')] = 0
    parameters['Gcphase'][('S', 'XY')] = 0
    parameters['Gcphase'][('S', 'XZ')] = 0
    
    parameters['Gcphase'][('S', 'YI')] = 5e-5*scale
    parameters['Gcphase'][('S', 'YX')] = 0
    parameters['Gcphase'][('S', 'YY')] = 0
    parameters['Gcphase'][('S', 'YZ')] = 0
    
    parameters['Gcphase'][('S', 'ZI')] = 5e-4*scale
    parameters['Gcphase'][('S', 'ZX')] = 0
    parameters['Gcphase'][('S', 'ZY')] = 0
    parameters['Gcphase'][('S', 'ZZ')] = 1e-3*scale
    
    # Parameter indexing, used throughout the analysis. It defines our "theta" vector
    return parameters

def make_s_error_param_dict(scale=1, include_loss=False):
    parameters = {}
    
    for gatename in ['Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2']:
        parameters[gatename] = {}
        parameters[gatename][('S', 'X')] = 0.2e-3*scale
        parameters[gatename][('S', 'Z')] = 0.2e-3*scale 
        parameters[gatename][('S', 'Y')] = 0.4e-3*scale

    parameters['Gi'] = {}
    parameters['Gi'][('S', 'X')] = 0
    parameters['Gi'][('S', 'Z')] = 0
    parameters['Gi'][('S', 'Y')] = 0
    
    parameters['Gcphase'] = {}
    parameters['Gcphase'][('S', 'IX')] = 5e-5*scale
    parameters['Gcphase'][('S', 'IY')] = 5e-5*scale
    parameters['Gcphase'][('S', 'IZ')] = 5e-4*scale
    
    parameters['Gcphase'][('S', 'XI')] = 5e-5*scale
    
    parameters['Gcphase'][('S', 'YI')] = 5e-5*scale

    parameters['Gcphase'][('S', 'ZI')] = (5e-4)*scale
    
    parameters['Gcphase'][('S', 'ZZ')] = (1e-3)*scale

    if include_loss:
        parameters['Gl'][('S','X')] = 5e-5*scale
        parameters['Gl'][('S','Y')] = 5e-5*scale
        parameters['Gl'][('S','Z')] = 1e-3*scale
    return parameters

def make_sh_sweep_error_param_dict(scale=1, h_param=0, spam_error=0):
    #make some portion of the error H error for some fixed error model
    s_parameters = {}
    h_parameters = {}
    
    for gatename in ['Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2']:
        s_parameters[gatename] = {}
        h_parameters[gatename] = {}
        s_parameters[gatename][('S', 'X')] = 0.1e-3*scale
        s_parameters[gatename][('S', 'Z')] = 0.4e-3*scale*(1-h_param)+0.1e-3*scale
        s_parameters[gatename][('S', 'Y')] = 0.1e-3*scale
        h_parameters[gatename][('H', 'Z')] = np.sqrt(0.4e-3*h_param*scale)


    s_parameters['Gi'] = {}
    s_parameters['Gi'][('S', 'X')] = 0
    s_parameters['Gi'][('S', 'Z')] = 0
    s_parameters['Gi'][('S', 'Y')] = 0

    h_parameters['Gi'] = {}
    h_parameters['Gi'][('H', 'X')] = 0
    h_parameters['Gi'][('H', 'Z')] = 0
    h_parameters['Gi'][('H', 'Y')] = 0

    s_parameters['Gc0'] = {}
    s_parameters['Gc0'][('S', 'X')] = 0.1e-3*scale
    s_parameters['Gc0'][('S', 'Z')] = 0.4e-3*scale*(1-h_param) + 0.1e-3*scale
    s_parameters['Gc0'][('S', 'Y')] = 0.1e-3*scale

    h_parameters['Gc0'] = {}
    h_parameters['Gc0'][('H', 'X')] = 0
    h_parameters['Gc0'][('H', 'Z')] = np.sqrt(0.4e-3*h_param*scale)
    h_parameters['Gc0'][('H', 'Y')] = 0
    
    s_parameters['Gcphase'] = {}
    h_parameters['Gcphase'] = {}
    s_parameters['Gcphase'][('S', 'IX')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'IY')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'IZ')] = 5e-4*scale*(1-h_param)+5e-5*scale
    
    s_parameters['Gcphase'][('S', 'XI')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'YI')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'ZI')] = (5e-4)*scale*(1-h_param)+5e-5*scale
    
    s_parameters['Gcphase'][('S', 'XX')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'XY')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'XZ')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'YX')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'YY')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'YZ')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'ZX')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'ZY')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'ZZ')] = (1e-3)*scale*(1-h_param)+5e-5*scale
    
    h_parameters['Gcphase'][('H', 'ZI')] = np.sqrt(5e-4*scale*h_param)
    h_parameters['Gcphase'][('H', 'IZ')] = np.sqrt(5e-4*scale*h_param)
    h_parameters['Gcphase'][('H', 'ZZ')] = np.sqrt((1e-3)*scale*h_param)

    s_parameters['Mdefault'] = {('S', 'X'): spam_error}
    s_parameters['rho0'] = {('S', 'X'): spam_error}

    return s_parameters, h_parameters

def make_sh_sweep_error_param_dict2(scale=1, h_param=0, spam_error=0):
    #make some portion of the error H error for some fixed error model
    s_parameters = {}
    h_parameters = {}
    
    for gatename in ['Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2']:
        s_parameters[gatename] = {}
        h_parameters[gatename] = {}
        s_parameters[gatename][('S', 'X')] = 0.1e-3*scale
        s_parameters[gatename][('S', 'Z')] = 0.4e-3*scale
        s_parameters[gatename][('S', 'Y')] = 0.1e-3*scale
        h_parameters[gatename][('H', 'Z')] = 0


    s_parameters['Gi'] = {}
    s_parameters['Gi'][('S', 'X')] = 0
    s_parameters['Gi'][('S', 'Z')] = 0
    s_parameters['Gi'][('S', 'Y')] = 0

    h_parameters['Gi'] = {}
    h_parameters['Gi'][('H', 'X')] = 0
    h_parameters['Gi'][('H', 'Z')] = 0
    h_parameters['Gi'][('H', 'Y')] = 0

    s_parameters['Gc0'] = {}
    s_parameters['Gc0'][('S', 'X')] = 0.1e-3*scale
    s_parameters['Gc0'][('S', 'Z')] = 0.4e-3*scale*(1-h_param) + 0.1e-3*scale
    s_parameters['Gc0'][('S', 'Y')] = 0.1e-3*scale

    h_parameters['Gc0'] = {}
    h_parameters['Gc0'][('H', 'X')] = 0
    h_parameters['Gc0'][('H', 'Z')] = np.sqrt(0.4e-3*h_param*scale)
    h_parameters['Gc0'][('H', 'Y')] = 0
    
    s_parameters['Gcphase'] = {}
    h_parameters['Gcphase'] = {}
    s_parameters['Gcphase'][('S', 'IX')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'IY')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'IZ')] = 5e-5*scale
    
    s_parameters['Gcphase'][('S', 'XI')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'YI')] = 5e-5*scale
    s_parameters['Gcphase'][('S', 'ZI')] = (5e-5)*scale
    
    s_parameters['Gcphase'][('S', 'XX')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'XY')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'XZ')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'YX')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'YY')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'YZ')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'ZX')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'ZY')] = (5e-5)*scale
    s_parameters['Gcphase'][('S', 'ZZ')] = (1e-3)*scale*(1-h_param)+5e-5*scale
    
    h_parameters['Gcphase'][('H', 'ZI')] = 0
    h_parameters['Gcphase'][('H', 'IZ')] = 0
    h_parameters['Gcphase'][('H', 'ZZ')] = np.sqrt((1e-3)*scale*h_param)

    s_parameters['Mdefault'] = {('S', 'X'): spam_error}
    s_parameters['rho0'] = {('S', 'X'): spam_error}

    return s_parameters, h_parameters