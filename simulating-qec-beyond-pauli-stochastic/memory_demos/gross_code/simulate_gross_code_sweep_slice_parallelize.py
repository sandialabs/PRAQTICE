#This script loads pre-computed DEMs svaed from running simulate_gross_code_sweep.py and generates more samples from them
import stim
import numpy as np
import pickle

from stimbposd import BPOSD, SinterDecoder_BPOSD, sinter_decoders

import time

from pygsti.extras.dem_construction import pygsti_object_builders as obj
from pygsti.extras.dem_construction import dem_tools as dems

import pygsti.tools.errgenproptools as eprop
from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator

from bb_tools import *
import multiprocessing
import itertools

import pauli_twirled_tools as pt
import pygsti

import os
import sys

module_path = os.path.abspath('../../../../BeamSearchDecoder')
sys.path.append(module_path)

import beamsearch

#load the relevant DEMs, repeatedly sample from them/decode, and save the results
h_scales1 = np.linspace(-1,1,11)
h_scales2 = [0.4, 0.6000000000000001]
rounds=2


shots = 100000

def simulate_bb(idx_and_hs):
    #load the relevant results
    #get the DEM
    i, (h1, h2) = idx_and_hs
    sim_name='_cancel'
    rounds=2
    
    with open(f'data/bb_code_simulation_{h1}_{h2}_{rounds}{sim_name}.pkl', 'rb') as f: 
        results_dict = pickle.load(f)
    f.close()

    dem=results_dict['dem']
    decoder = beamsearch.BeamSearch(dem)
    
    sampler = dem.compile_sampler()
    samples = sampler.sample(shots=shots)
    
    predicted_observables = decoder.decode_batch(samples[0])
    num_mistakes = np.sum(np.any(predicted_observables != samples[1], axis=1))
    
    #get the pt_dem
    pt_dem=results_dict['pt_dem']
    
    pt_decoder = beamsearch.BeamSearch(pt_dem)
    
    pt_predicted_observables = pt_decoder.decode_batch(samples[0])
    pt_num_mistakes = np.sum(np.any(pt_predicted_observables != samples[1], axis=1))

    pt_sampler = pt_dem.compile_sampler()
    pt_samples = pt_sampler.sample(shots=shots)

    pt_pt_predicted_observables = pt_decoder.decode_batch(pt_samples[0])
    pt_pt_num_mistakes = np.sum(np.any(pt_pt_predicted_observables != pt_samples[1], axis=1))
    

    print(i, h1, h2, num_mistakes, pt_num_mistakes, pt_pt_num_mistakes)
    
    ###SAVE RESULTS###
    results_dict = {'samples': samples,
                    'obs': predicted_observables,
                    'pt_obs': pt_predicted_observables,
                    'pt_samples': pt_samples, 
                    'pt_pt_obs': pt_pt_predicted_observables}

#save samples, predicted observables in a relevant file
    
    with open(f'data/bb_code_simulation_{h1}_{h2}_{rounds}_cancel_slice.pkl', 'wb+') as f:
        pickle.dump(results_dict, f)
    f.close()


if __name__=="__main__":
    for i in range(4,14):
        print(f'starting round {i}')
        with multiprocessing.Pool(processes=10) as pool:
                results = pool.map(simulate_bb, list((i,hs) for hs in itertools.product(h_scales1,h_scales2)))
