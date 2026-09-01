###goal: calculate LER for surface code with coherent error###s
###sweep relative contribution of S and H error and look at change in logical error rate. 
###also compare Pauli-twirl informed decoder and the approximate DEM decoder
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

from pygsti.baseobjs.statespace import QubitSpace as _QubitSpace
import pygsti.baseobjs as _bo
from pygsti.baseobjs.errorgenlabel import GlobalElementaryErrorgenLabel as _GlobalElementaryErrorgenLabel, \
                                          LocalElementaryErrorgenLabel as _LocalElementaryErrorgenLabel
import pygsti.tools.lindbladtools as lbd

from itertools import combinations, product

import pymatching
ds=[3,5,7] 
ps = np.linspace(1.8,10.8, num=11) 
n_shots = 10000000

num_shca_terms = [10,3,2,2] #number of terms of each EEG type in the models for 2Q gates. In the order S, H, C, A
s_params=np.array([0.001/15, 0.0002/15])#parameters for sampling 2Q gate errors
h_params=np.array([0,np.sqrt(0.001)/5])

spam_error = 0.001
num_models = 20
from stim_circuits_for_pt_sims import *

def generate_random_sparse_shca_error_model(nqs, num_shca_terms, s_params, h_params, seed=1234):
    rng = np.random.default_rng(seed)
    state_space = _QubitSpace.cast(nqs)
    
    #create an error generator basis according the our weight specs
    errorgen_basis = _bo.CompleteElementaryErrorgenBasis('PP', state_space, elementary_errorgen_types=['S','H'],
                                                         default_label_type='local')
    
    #Get the labels, broken out by sector, of each of the error generators in this basis.
    errgen_labels_H = lbd._sort_errorgen_labels(errorgen_basis.sublabels('H'))
    errgen_labels_S = lbd._sort_errorgen_labels(errorgen_basis.sublabels('S'))
    
    edict = {}
    sa_edict = {}
    sc_edict = {}
    #pick h error gens and sample rates
    chosen_labels_H = rng.choice(errgen_labels_H, num_shca_terms[1] ,replace=False)
    num_H_rates = num_shca_terms[1]
    edict.update({lbl: sign*val for lbl,sign, val in zip(errgen_labels_H, rng.choice([-1, 1], size = num_H_rates), rng.normal(loc=h_params[0], scale=h_params[1], size = num_H_rates))})
    #pick s error gens and generate random s rates
    num_S_rates = num_shca_terms[0]
    chosen_labels_S = [errgen_labels_S[i] for i in rng.choice(len(errgen_labels_S), num_shca_terms[0] ,replace=False)]  
    random_S_vals = [np.abs(j) for j in rng.normal(loc=s_params[0], scale=s_params[1], size = num_S_rates)]

    random_S_vals1 = [np.abs(j) for j in rng.normal(loc=s_params[0], scale=s_params[1], size = num_S_rates)]
    random_S_vals2 = [np.abs(j) for j in rng.normal(loc=s_params[0], scale=s_params[1], size = num_S_rates)]
    edict.update({lbl: val1+val2 for lbl,val1,val2 in zip(chosen_labels_S, random_S_vals1, random_S_vals2)})
    
    sa_edict.update({lbl: val for lbl,val in zip(chosen_labels_S, random_S_vals1)})
    sc_edict.update({lbl: val for lbl,val in zip(chosen_labels_S, random_S_vals2)})

    #I just want to generate the valid C and A labels, given some valid S labels
    allowed_paulis = [eeg.basis_element_labels[0] for eeg in chosen_labels_S]
    pauli_pairs = list(combinations(allowed_paulis, 2))
    
    chosen_paulis_C = [pauli_pairs[i] for i in rng.choice(len(pauli_pairs), num_shca_terms[2], replace=False)]
    
    chosen_labels_C = []
    for ps in chosen_paulis_C:
        lbl = _LocalElementaryErrorgenLabel('C', ps)
        chosen_labels_C.append(lbl)
        #choose their rates
        max_rate = np.sqrt(sc_edict[_LocalElementaryErrorgenLabel('S', (ps[0],))]*sc_edict[_LocalElementaryErrorgenLabel('S', (ps[1],))])
        rate = rng.choice([-1, 1])*np.random.uniform(0, max_rate)
        edict[lbl] = rate
    
    chosen_paulis_A = [pauli_pairs[i] for i in rng.choice(len(pauli_pairs), num_shca_terms[3], replace=False)]
    
    chosen_labels_A = []
    for ps in chosen_paulis_A:
        lbl = _LocalElementaryErrorgenLabel('A', ps)
        chosen_labels_A.append(lbl)
        #choose their rates
        max_rate = np.sqrt(sa_edict[_LocalElementaryErrorgenLabel('S', (ps[0],))]*sa_edict[_LocalElementaryErrorgenLabel('S', (ps[1],))])
        rate = rng.choice([-1, 1])*np.random.uniform(0, max_rate)
        edict[lbl] = rate
    
    return edict

def simulate_random_err(idx):
    lers = {d:[] for d in ds}
    s_lers = {d:[] for d in ds}

    error_rates_dict = {'Gcnot': generate_random_sparse_shca_error_model(2, num_shca_terms, s_params, h_params, seed=idx), 
           'Gh': {('S', 'X'): s_params[0]*0.1/3,('S', 'Y'): s_params[0]*0.1/3,('S', 'Z'): s_params[0]*0.1/3}}
    error_rates_dict['Mdefault'] = {('S', 'X'): spam_error}
    error_rates_dict['rho0'] = {('S', 'X'): spam_error}
    
    for d in ds:
        for p in ps:
            rounds = d
            h_error_rates_dict = {}
            for g,edict in error_rates_dict.items():
                #print(edict)
                if g=='Gcnot':
                    h_error_rates_dict[g] = {k:v*p if k.errorgen_type!='H' else v*np.sqrt(p) for k,v in edict.items()}
                else:
                    h_error_rates_dict[g] = {k:v*p if k[0]!='H' else v*np.sqrt(p) for k,v in edict.items()}
            
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
    
            h_model = obj.build_model(h_error_rates_dict, mr_pspec, oneQ_gate_names=['Gh'], twoQ_gate_names=['Gcnot'])
            
            time1 = time.time()
            h_egp = ErrorGeneratorPropagator(h_model)
            eoc_eeg = h_egp.propagate_errorgens_bch(mr_pcircuit, bch_order=1, include_spam=False)
            time2 = time.time()
            print(f'{time2-time1} to propagate errors for d={d}')

            total_dem = dems.generate_dem_higher_order(dets_as_pauli_strings, eoc_eeg, sim, zassenhaus_order=1, add_type='add')
    
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
        
    
            with open(f'data/surface_code_{d}_random_models_{idx}_strength_{p}.pkl', 'wb+') as f:
                pickle.dump([ler[0],total_dem, h_error_rates_dict], f)
            f.close()

    print(f'completed d={d}')

    return lers




if __name__=="__main__":
    with multiprocessing.Pool(processes=5) as pool:
        results = pool.map(simulate_random_err, range(num_models))
