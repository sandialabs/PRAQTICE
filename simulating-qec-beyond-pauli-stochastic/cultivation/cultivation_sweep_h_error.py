'''
simulations sweeping the %h error in a simple error model and looking at (1) the logical error rate, and (2) the discard rate

currently coded for d=3 first two stages of cultivation. 
improvements will 
'''

import numpy as np
import stim
import re
import pygsti
import time
import itertools
import multiprocessing
import pickle

from collections import Counter

import pygsti.tools.errgenproptools as eprop
from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator

from cultiv_circuits_cleaned import *
from error_models import *


#this code borrows a few functions and circuits from Gidney et al.'s code for cultivation
import pathlib
import sys
src_path = pathlib.Path('../cultivation/cultivation-paper-code/code/src')
assert src_path.exists()
sys.path.append(str(src_path))

import cultiv
import gen

#As of 04-26-26, this code requires being on the pyGSTi branch feature-dems-from-errorgens
from pygsti.extras.dem_construction import pygsti_object_builders as obj
from pygsti.extras.dem_construction import dem_tools as dems

num_processes = 20
n_shots = 100000000 #set this number much larger to get good statistics for relevant error rates; small for testing
n_runs = 100 #run more times to collect more samples
h_fractions = np.linspace(1,0, num=20) #np.linspace(1,0, num=20) #fraction of generator infidelity allocated to coherent error 
strength = 1#1.6#1.4 #overall scale factor for gate generator infidelities in the error model

incl_spam = True
#rate of bit flip error on each qubit for prep/measure
if incl_spam:
    spam_error = 0.001 #*strength#0.0015
else:
    spam_error = 0 

###########################
def create_processor_spec(pcircuit, qubit_labels):
    """
    Creates a processor spec that from a pyGSTi syndrome extraction circuit, containing the 
    required gates (CNOTs between qubits that are coupled in those circuits). Will then be
    used to define the noise model we simulate.
    """
    pyg_cstr = pcircuit.str
    # Find all the CNOT gates used in the circuit, using brute force string comprehension
    cleaned_pcstr = pyg_cstr.replace(']', '').replace('[', '').split('@')[0]
    connections = [(int(s.split(':')[1]), int(s.split(':')[2])) for s in cleaned_pcstr.split('G') if len(s.split(':')) > 2]

    availability={'Gcphase':connections, 'Gcnot':connections}
    pspec = pygsti.processors.QubitProcessorSpec(len(qubit_labels), ['Gcphase', 'Gcnot' ,'Gh', 'Gypi2', 'Gympi2', 'Gxpi2', 'Gxpi', 'Gxmpi2', 'Gi','Gc0'], qubit_labels=qubit_labels,
                                             availability=availability)
    
    return pspec

#######BASIC SETUP#########
#inject_style should be 'unitary' or 'degenerate' for these simulations. 

def run_simulation(h_fraction):
    for inject_style in ['unitary','degenerate']:
        if inject_style=='unitary':
            cultiv_c_basic = create_d3_cultiv_circ()
            cultiv_circuit = gen.transpile_to_z_basis_interaction_circuit(cultiv_c_basic) 
        
        else:
            cultiv_c_basic =  cultiv.make_inject_and_cultivate_circuit(inject_style=inject_style, basis='X', dcolor=3)
            cultiv_circuit = gen.transpile_to_z_basis_interaction_circuit(cultiv_c_basic) 

        c, qubit_mapping, measurements, detectors = obj.stim_to_pygsti_circuit(cultiv_circuit, range(cultiv_circuit.num_qubits), qubit_relabelling_dict=None, show_qubit_mappings=False, include_idles=True, include_observables=True)
        
        num_qubits = len(c.line_labels)
        qubit_labels = range(num_qubits)
        
        pspec = create_processor_spec(c, qubit_labels)
        
        twoq_gates = ['Gcphase']
        oneq_gates = [g for g in pspec.gate_names if g != 'Gcphase' and g != 'Gcnot']
        
        s_error_rates_dict, h_error_rates_dict = make_sh_sweep_error_param_dict(scale=strength, h_param=h_fraction, spam_error=spam_error)
    
        h_model = obj.build_model(h_error_rates_dict, pspec, oneq_gates, twoq_gates, meas_err=0, prep_err=0)
        h_egp = ErrorGeneratorPropagator(h_model)
    
        s_model = obj.build_model(s_error_rates_dict, pspec, oneq_gates, twoq_gates, meas_err=spam_error, prep_err=spam_error)
        s_egp = ErrorGeneratorPropagator(s_model)
    
        serial_c = c.serialize()

    
        stim_c_no_extras = obj.pygsti_c_to_stim(serial_c)
    
        tableau = stim.Tableau.from_circuit(stim_c_no_extras, ignore_measurement=True) 
        inverse_tableau = tableau.inverse()
    
        dps = [dems.get_detector_as_parity(d, measurements, num_qubits) for d in detectors]
    
        anc_qubits = [m[1][0] for m in measurements if m[0]=='Z']

        sim = stim.TableauSimulator()
        sim.set_inverse_tableau(inverse_tableau)
        
        s_egp = ErrorGeneratorPropagator(s_model)
        h_egp = ErrorGeneratorPropagator(h_model)
        time1 = time.time()
        eoc_s_eeg = s_egp.propagate_errorgens_bch(serial_c, bch_order=1, include_spam=incl_spam)
        eoc_h_eeg = h_egp.propagate_errorgens_bch(serial_c, bch_order=1, include_spam=False)
        time2 = time.time()
        print(f'{time2-time1} to propagate S errors')
    
        h_dem = dems.generate_dem(dps, eoc_h_eeg, sim, anc_qubits)
        s_dem = dems.generate_dem(dps, eoc_s_eeg, sim, anc_qubits)
    
        dem = Counter(h_dem)+Counter(s_dem)
    
        print(f'generated DEM for {h_fraction}')
    
        stim_dem_str = dems.format_dem_stim(dem, n_logical=0)
        stim_dem = stim.DetectorErrorModel(stim_dem_str)
        h_sampler = stim_dem.compile_sampler()
        n_dets = len(list(dem.keys())[0])

        
        lers = []
        p_keeps = []
        for j in range(n_runs):
            print(f'on run {j}')
            h_samples = h_sampler.sample(shots=n_shots)[0]
            h_dem_stats = Counter([''.join(['1' if b else '0' for b in row]) for row in h_samples])
    
            print(f'sampled DEM for {h_fraction}')
    
            ler = h_dem_stats['0'*(n_dets-1)+'1']/n_shots
    
            p_keep = (h_dem_stats['0'*(n_dets-1)+'1']+h_dem_stats['0'*n_dets])/n_shots
            lers.append(ler)
            p_keeps.append(p_keep)

            print(h_fraction,ler,p_keep)
        
        # ler = sum(lers)/n_runs
        # p_keep = sum(p_keeps)/n_runs
        
            #record results
            results = {'frac_kept': p_keep, 'ler': ler, 'dem': dem, 'h_frac':h_fraction,'strength':strength, 'spam_error':spam_error, 'error_rates_dicts':(s_error_rates_dict, h_error_rates_dict), 'n_runs':j, 'lers': lers, 'pkeeps':p_keep, 'h_dem':h_dem, 's_dem':Counter(s_dem)}
            
            with open(f'cultivation_data/cultiv_h_sweep_{inject_style}_{h_fraction}_{strength}_run_{j+50}_test.pkl', 'wb+') as f:
                pickle.dump(results, f)
            f.close()

        print(f'saved results for {h_fraction} {inject_style}')
    
    return None

if __name__=="__main__":
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(run_simulation, h_fractions)
    
