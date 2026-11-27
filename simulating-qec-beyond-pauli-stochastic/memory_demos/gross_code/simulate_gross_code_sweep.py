import stim
import numpy as np
import pickle
import pygsti

import time

#As of 04-26-26, this code requires being on the pyGSTi branch feature-dems-from-errorgens
from pygsti.extras.dem_construction import pygsti_object_builders as obj
from pygsti.extras.dem_construction import dem_tools as dems

import pygsti.tools.errgenproptools as eprop
from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator

from bb_tools import *
import multiprocessing
import itertools

import pauli_twirled_tools as pt

#running this requires the Beam Search Decoder, but this can be repleased with BPOSD or a similar BP-like decoder.
import os
import sys
module_path = os.path.abspath('../../BeamSearchDecoder')
sys.path.append(module_path)

import beamsearch

h_scales1 = np.linspace(-1,1,11)
h_scales2 = np.linspace(0,1,11)
rounds=2

def make_sh_sweep_error_param_dict(p=1, h_scale=0, h_scale_idle=0, spam_error=0):
    h_error_rates_dict = {}
    h_error_rates_dict['Gi'] = {('H','X'):np.sqrt(0.001*p*h_scale_idle)}
    h_error_rates_dict['Gcnot'] = {('H','IX'): np.sqrt(0.001*p*h_scale),
                                ('H', 'ZI'):np.sqrt(0.001*p*h_scale),
                                ('H', 'ZX'):np.sqrt(0.001*p*h_scale)}
    h_error_rates_dict['Gi'].update({('S','X'):0.0005*p+0.001*p*(1-h_scale_idle),('S','Y'):0.0005*p,('S','Z'):0.0005*p})
    h_error_rates_dict['Gh'] = {('S','Z'):0.001*p,('S','Y'):0.001*p,('S','X'):0.001*p}
    h_error_rates_dict['Gcnot'].update({('S','IX'): 0.001*p/15+0.001*p*(1-h_scale),('S', 'IY'):0.001*p/15,('S', 'IZ'):0.001*p/15,
                                  ('S', 'XX'):0.001*p/15,('S', 'XY'):0.001*p/15,('S', 'XZ'):0.001*p/15,
                                   ('S', 'YX'):0.001*p/15, ('S', 'YY'):0.001*p/15, ('S', 'YZ'):0.001*p/15,
                                  ('S', 'ZX'):0.001*p/15+0.001*p*(1-h_scale), ('S', 'ZY'):0.001*p/15, ('S', 'ZZ'):0.001*p/15, 
                                       ('S','XI'): 0.001*p/15,('S', 'YI'):0.001*p/15,('S', 'ZI'):0.001*p/15+0.001*p*(1-h_scale)})

    h_error_rates_dict['Mdefault'] = {('S', 'X'): spam_error}
    h_error_rates_dict['rho0'] = {('S', 'X'): 0}

    return h_error_rates_dict

def make_sh_cancel_sweep_error_param_dict2(p=1, h_scale=0, h_scale_idle=0, spam_error=0):
    h_error_rates_dict = {}
    h_error_rates_dict['Gi'] = {('S','X'):0.001*p,('S','Y'):0.001*p,('S','Z'):0.001*p, ('H','X'):np.sqrt(0.0003*p)*h_scale_idle}
    h_error_rates_dict['Gcnot'] = {('H', 'XX'):np.sqrt(0.001*p)*h_scale}
    h_error_rates_dict['Gh'] = {('S','Z'):0.001*p,('S','Y'):0.001*p,('S','X'):0.001*p}
    h_error_rates_dict['Gcnot'].update({('S','IX'): 0.001*p,('S', 'IY'):0.001*p,('S', 'IZ'):0.001*p,
                                       ('S','XI'): 0.001*p,('S', 'YI'):0.001*p,('S', 'ZI'):0.001*p})

    h_error_rates_dict['Mdefault'] = {('S', 'X'): spam_error}
    h_error_rates_dict['rho0'] = {('S', 'X'): 0}


    return h_error_rates_dict

circuit = build_gross_bb_repeated(rounds=rounds)

c, qubit_mapping, measurements, detectors = obj.stim_to_pygsti_circuit(circuit, range(circuit.num_qubits), qubit_relabelling_dict=None, show_qubit_mappings=False, include_idles=False, include_meas_idles=True, include_observables=True)
mr_pspec = obj.create_processor_spec(c, c.line_labels, gates=['Gh','Gi','Gcnot']) 
n_qubits = len(c.line_labels)

serial_c = c.serialize()
stim_c_no_extras = obj.pygsti_c_to_stim(serial_c)

tableau = stim.Tableau.from_circuit(stim_c_no_extras, ignore_measurement=True) 
inverse_tableau = tableau.inverse()

sim = stim.TableauSimulator()
sim.set_inverse_tableau(inverse_tableau)

spam_error=0.001

shots = 1000000

def simulate_bb(hs):
    h1, h2 = hs
    #build error model
    h_error_rates_dict = make_sh_cancel_sweep_error_param_dict2(p=1, h_scale=h1, h_scale_idle=h2, spam_error=spam_error) 
    h_model = obj.build_model(h_error_rates_dict, mr_pspec, oneQ_gate_names=['Gh','Gi'], twoQ_gate_names=['Gcnot'])
    
    dets_as_pauli_strings = [dems.get_detector_as_parity(d, measurements, n_qubits) for d in detectors]
    
    ###Errorprop setup
    h_egp = ErrorGeneratorPropagator(h_model)
    time1 = time.time()
    eoc_eeg = h_egp.propagate_errorgens_bch(serial_c, bch_order=1, include_spam=True)
    time2 = time.time()
    print(f'{time2-time1} to propagate errors')
    
    total_dem = dems.generate_dem_higher_order(dets_as_pauli_strings, eoc_eeg, sim, zassenhaus_order=1, add_type='add')
    dem_str = dems.format_dem_stim(total_dem, n_logical=12)
    
    dem = stim.DetectorErrorModel(dem_str)
    decoder = beamsearch.BeamSearch(dem)
    
    sampler = dem.compile_sampler()
    samples = sampler.sample(shots=shots)
    
    predicted_observables = decoder.decode_batch(samples[0])
    num_mistakes = np.sum(np.any(predicted_observables != samples[1], axis=1))
    
    #####Pauli Twirled#####
    error_pm_dict = {'Gcnot': h_model.operation_blks['gates'][pygsti.baseobjs.label.Label('Gcnot',(147,216))].factorops[1].to_dense(),
                    'Gh': h_model.operation_blks['gates'][pygsti.baseobjs.label.Label('Gh',(147))].factorops[1].to_dense(),
                    'Gi': h_model.operation_blks['gates'][pygsti.baseobjs.label.Label('Gi',(71))].factorops[1].to_dense()}
    
    pt_dem = pt.build_pauli_twirled_model(error_pm_dict, circuit, spam_error=spam_error, include_idles=False, include_meas_idles=True)
    
    pt_decoder = beamsearch.BeamSearch(pt_dem)
    
    pt_predicted_observables = pt_decoder.decode_batch(samples[0])
    pt_num_mistakes = np.sum(np.any(pt_predicted_observables != samples[1], axis=1))

    pt_sampler = pt_dem.compile_sampler()
    pt_samples = pt_sampler.sample(shots=shots)

    pt_pt_predicted_observables = pt_decoder.decode_batch(pt_samples[0])
    pt_pt_num_mistakes = np.sum(np.any(pt_pt_predicted_observables != pt_samples[1], axis=1))
    

    print(h1, h2, num_mistakes, pt_num_mistakes, pt_pt_num_mistakes)
    
    ###SAVE RESULTS###
    results_dict = {'dem': dem,
                    'pt_dem': pt_dem,
                    'samples': samples,
                    'obs': predicted_observables,
                    'pt_obs': pt_predicted_observables,
                    'pt_samples': pt_samples, 
                    'pt_pt_obs': pt_pt_predicted_observables}
    
    with open(f'../data/bb_code_simulation_{h1}_{h2}_{rounds}_cancel.pkl', 'wb+') as f:
        pickle.dump(results_dict, f)
    f.close()


if __name__=="__main__":
    with multiprocessing.Pool(processes=10) as pool:
        results = pool.map(simulate_bb, list(itertools.product(h_scales1,h_scales2)))
