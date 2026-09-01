import numpy as np
import stim
import pygsti
import time
import itertools

import multiprocessing
from collections import Counter
import pickle

import pygsti.tools.errgenproptools as eprop
from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator

#As of 04-26-26, this code requires being on the pyGSTi branch feature-dems-from-errorgens
from pygsti.extras.dem_construction import pygsti_object_builders as obj
from pygsti.extras.dem_construction import dem_tools as dems

from qec_errorprop.codes import surface_code_setup as surface

import pauli_twirled_tools as pt 

import pymatching
ds=[3,5,7]
ps = np.linspace(0.3,0.5, num=5) 
ps = np.append(ps, np.linspace(0.006,0.3,num=15))
ps = np.append(ps, np.linspace(0.1,0.3,num=8))
hscales = np.linspace(0,1, num=5)
n_shots = 10000000

from stim_circuits_for_pt_sims import *

def make_sh_sweep_error_param_dict(p=1, h_scale=0, spam_error=0):
    h_error_rates_dict = {}
    h_error_rates_dict['Gh'] = {('H','Z'):np.sqrt(0.0001*p*h_scale),('H','X'):np.sqrt(0.0001*p*h_scale)}
    h_error_rates_dict['Gcnot'] = {('H','IX'): np.sqrt(0.01*p*h_scale),
                                ('H', 'ZI'):np.sqrt(0.01*p*h_scale),
                                ('H', 'ZX'):np.sqrt(0.01*p*h_scale)}
    s_error_rates_dict = {}
    s_error_rates_dict['Gh'] = {('S','Z'):0.0001*p+0.0001*p*(1-h_scale),('S','Y'):0.0001*p,('S','X'):0.0001*p*(1-h_scale)}
    s_error_rates_dict['Gcnot'] = {('S','IX'): 0.01*p/15+0.01*p*(1-h_scale),('S', 'IY'):0.01*p/15,('S', 'IZ'):0.01*p/15,
                                  ('S', 'XX'):0.01*p/15,('S', 'XY'):0.01*p/15,('S', 'XZ'):0.01*p/15,
                                   ('S', 'YX'):0.01*p/15, ('S', 'YY'):0.01*p/15, ('S', 'YZ'):0.01*p/15,
                                  ('S', 'ZX'):0.01*p/15+0.01*p*(1-h_scale), ('S', 'ZY'):0.01*p/15, ('S', 'ZZ'):0.01*p/15, 
                                  ('S','XI'): 0.01*p/15,('S', 'YI'):0.01*p/15,('S', 'ZI'):0.01*p/15+0.01*p*(1-h_scale)}

    s_error_rates_dict['Mdefault'] = {('S', 'X'): spam_error}
    s_error_rates_dict['rho0'] = {('S', 'X'): spam_error}

    return s_error_rates_dict, h_error_rates_dict

def simulate_sh_err(p):
    lers = {d:[] for d in ds}
    s_lers = {d:[] for d in ds}
    
    for d in ds:
        for hscale in hscales:
            rounds = d
            s_error_rates_dict, h_error_rates_dict = make_sh_sweep_error_param_dict(p=p, h_scale=hscale, spam_error=0)
            
            stimc = stim.Circuit.generated("surface_code:rotated_memory_z", distance=d, rounds=rounds)
            stimc = pt.unroll_repeat_blocks(stimc)
            mr_pcircuit, qubit_mapping, measurements, detectors = obj.stim_to_pygsti_circuit(stimc, range(stimc.num_qubits)
                                                                                             , qubit_relabelling_dict=None, show_qubit_mappings=False, include_idles=False, include_meas_idles=False, include_observables=True)
            mr_pspec = obj.create_processor_spec(mr_pcircuit, mr_pcircuit.line_labels, gates=['Gh','Gcnot']) 

            n_qubits = len(mr_pcircuit.line_labels)
            
            dets_as_pauli_strings = [dems.get_detector_as_parity(d, measurements, n_qubits) for d in detectors]

            serial_c = mr_pcircuit.serialize()
            stim_c_no_extras = obj.pygsti_c_to_stim(serial_c)
            tableau = stim.Tableau.from_circuit(stim_c_no_extras, ignore_measurement=True) 
            inverse_tableau = tableau.inverse()
            
            sim = stim.TableauSimulator()
            sim.set_inverse_tableau(inverse_tableau)

            time1 = time.time()
            
            mr_pcircuit = mr_pcircuit.serialize()
            for k in h_error_rates_dict:
                h_error_rates_dict[k].update(s_error_rates_dict[k])
    
            h_model = obj.build_model(h_error_rates_dict, mr_pspec, oneQ_gate_names=['Gh'], twoQ_gate_names=['Gcnot'])

            time1 = time.time()
            h_egp = ErrorGeneratorPropagator(h_model)
            eoc_eeg = h_egp.propagate_errorgens_bch(mr_pcircuit, bch_order=1, include_spam=False)
            time2 = time.time()
            print(f'{time2-time1} to propagate errors for d={d}')

            total_dem = dems.generate_dem_higher_order(dets_as_pauli_strings, eoc_eeg, sim, zassenhaus_order=1, add_type='add')
        
            with open(f'data/sh_surface_code_dem_{d}_{p}.pkl', 'wb+') as f:
                pickle.dump(total_dem,f)
            f.close()
    
            print(f'created DEM for d={d}')

            stim_dem_str = dems.format_dem_stim(total_dem, n_logical=1)
            stim_dem = stim.DetectorErrorModel(stim_dem_str)

            dec = pymatching.Matching.from_detector_error_model(stim_dem)
            
            h_sampler = stim_dem.compile_sampler()
            h_samples = h_sampler.sample(shots=n_shots)[0]
    
            syndrome, log_flips, _ = h_sampler.sample(shots=n_shots)
            
            predictions = dec.decode_batch(syndrome)
            log_errors = (predictions != log_flips)
    
            ler = sum(log_errors)/n_shots
            print(f'LER of {ler[0]} for d={d}, p={p}')
            lers[d].append(ler[0])
        

            ########PAULI-TWIRLED MODEL'S DEM################

            error_pm_dict = {'Gcnot': h_model.operation_blks['gates'][pygsti.baseobjs.label.Label('Gcnot',(2,3))].factorops[1].to_dense(),
                'Gh': h_model.operation_blks['gates'][pygsti.baseobjs.label.Label('Gh',(2))].factorops[1].to_dense(),
                'Gi': np.eye(4)}

            pt_c = stim.Circuit()
            for i in range(rounds):
                pt_c.append(circuit_to_add_noise_base[d])
            pt_c.append(eoc_meas[d])
            s_only_dem = pt.build_pauli_twirled_model(error_pm_dict, stimc)#the pauli twirl model DEM
            
            sdec = pymatching.Matching.from_detector_error_model(s_only_dem)
            s_predictions = sdec.decode_batch(syndrome)
            s_log_errors = (s_predictions != log_flips)
    
            s_ler = sum(s_log_errors)/n_shots
            s_lers[d].append(s_ler[0])
            print(f'LER of {s_ler[0]} for d={d}, p={p} with S decoder')
    
            with open(f'data/sh_surface_code_{d}_lers_h_sweep_{p}_{hscale}.pkl', 'wb+') as f:
                pickle.dump([ler[0],s_ler[0],total_dem, s_only_dem, h_error_rates_dict], f)
            f.close()

    print(f'completed d={d}')

    return lers




if __name__=="__main__":
    with multiprocessing.Pool(processes=10) as pool:
        results = pool.map(simulate_sh_err, ps)
