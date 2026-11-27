#this script creates DEMs at first and second order in approximations for surface code SE with random H error models#
#the goal is to compare different approximations#
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['GOTO_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1' 
os.environ['NUMEXPR_NUM_THREADS'] = '1' 
os.environ['VECLIB_MAXIMUM_THREADS'] = '1' 
os.environ['MKL_NUM_THREADS'] = '1' 

from mpi4py import MPI 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from pygsti.modelmembers.operations import ExpErrorgenOp
from pygsti.baseobjs.statespace import QubitSpace
import pauli_twirled_tools as pt
import stim
import pygsti
import pickle
import numpy as np
import scipy as sp
n_models = 50


d=3
rounds=2


def compute_outcome_distribution_from_dem(stim_dem):
    """
    Compute the outcome distribution from a DEM using log and Hadamard transform.

    Parameters:
    - dem: a detector error model

    Returns: 
    - prob_estimate: array of 2^n probabilities, in increasing binary order
    """
    # Convert DEM to attenuations
    attenuations = np.zeros(2**stim_dem.num_detectors, dtype=float)
    all_bitstrings = np.array([[int(bit) for bit in format(n, f'0{stim_dem.num_detectors}b')]for n in range(2**stim_dem.num_detectors)])
    for event in stim_dem:
        #print(event.args_copy())
        if len(event.args_copy())==1: #hack
            prob = event.args_copy()[0]
            targets = [target.val for target in event.targets_copy()]
            event = [1 if (stim_dem.num_detectors-idx-1) in targets else 0 for idx in range(stim_dem.num_detectors)]    
            if not 1-2*prob>0: print(prob, event, event.args_copy())
            attenuation = -np.log(1-2*prob)
            attenuations += attenuation * (1-(-1)**np.dot(all_bitstrings, event))/2
    
    # Compute polarizations from attenuations
    polarizations = np.exp(-attenuations)
    polarizations[0] = 1

    # Compute probabilities from polarizations
    probabilities = sp.linalg.hadamard(2**stim_dem.num_detectors) @ polarizations / 2**stim_dem.num_detectors
    nds = stim_dem.num_detectors

    probs_dict = {str(bin(b)[2:]).zfill(nds)[::-1]:p for b,p in enumerate(probabilities)}
    return probs_dict

def compute_probs_from_dem(dem):
    '''
    generates samples from a detector error model
    dem (dict): A dictionary specifying the DEM
    n_shots (int): number of samples to take
    '''
    stim_dem_str = dems.format_dem_stim(dem, n_logical=0)
    stim_dem = stim.DetectorErrorModel(stim_dem_str)
    est_probs = compute_all_probs(stim_dem)
    
    return est_probs

init_perfect_circuit = stim.Circuit('''R 1 3 5 8 10 12 15 17 19 2 9 11 13 14 16 18 25
    TICK
    H 2 11 16 25
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18
    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18
    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13
    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13
    TICK
    H 2 11 16 25
    TICK
    MR 2 9 11 13 14 16 18 25''')

circuit_to_add_noise = stim.Circuit('''H 2 11 16 25
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18
    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18
    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13
    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13
    TICK
    H 2 11 16 25
    TICK
    MR 2 9 11 13 14 16 18 25
    DETECTOR(2, 0, 0) rec[-8] rec[-16]
    DETECTOR(2, 2, 0) rec[-7] rec[-15]
    DETECTOR(4, 2, 0) rec[-6] rec[-14]
    DETECTOR(6, 2, 0) rec[-5] rec[-13]
    DETECTOR(0, 4, 0) rec[-4] rec[-12]
    DETECTOR(2, 4, 0) rec[-3] rec[-11]
    DETECTOR(4, 4, 0) rec[-2] rec[-10]
    DETECTOR(4, 6, 0) rec[-1] rec[-9]
    TICK
    H 2 11 16 25
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18
    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18
    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13
    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13
    TICK
    H 2 11 16 25
    TICK
    MR 2 9 11 13 14 16 18 25
    DETECTOR(2, 0, 0) rec[-8] rec[-16]
    DETECTOR(2, 2, 0) rec[-7] rec[-15]
    DETECTOR(4, 2, 0) rec[-6] rec[-14]
    DETECTOR(6, 2, 0) rec[-5] rec[-13]
    DETECTOR(0, 4, 0) rec[-4] rec[-12]
    DETECTOR(2, 4, 0) rec[-3] rec[-11]
    DETECTOR(4, 4, 0) rec[-2] rec[-10]
    DETECTOR(4, 6, 0) rec[-1] rec[-9]
    M 1 3 5 8 10 12 15 17 19
    DETECTOR rec[-7] rec[-8] rec[-9]''') #OBSERVABLE_INCLUDE(0) rec[-7] rec[-8] rec[-9]

def surface_code_create_pt_dem(error_dict):
    #process the dictionary of error rates into required format
    edict_processed = {}
    for gname, gdict in error_dict.items():
        edict_processed[gname] = {}
        for k,v in gdict.items():
            edict_processed[gname][(k[0],str(k[1]))] = v

    pm_edict = {'Gh':ExpErrorgenOp(pygsti.modelmembers.operations.LindbladErrorgen.from_elementary_errorgens(edict_processed['Gh'], state_space=QubitSpace(1))).to_dense(),
            'Gcnot':ExpErrorgenOp(pygsti.modelmembers.operations.LindbladErrorgen.from_elementary_errorgens(edict_processed['Gcnot'], state_space=QubitSpace(2))).to_dense(),
            'Gi':np.eye(4)}

    pt_dem = pt.build_pauli_twirled_model(pm_edict, circuit_to_add_noise, init_circuit=init_perfect_circuit)

    return pt_dem

n_models=5
erates = np.linspace(0.002, 0.02, num=11)
erates = np.append(erates,[0.001,0.003])

for erate in erates:
    for i in range(n_models):
        print(f'{erate} model {i}')
        try:
            with open(f'data/sweep_rounds_{rounds}_infidelity_{erate}_{i}.txt', 'rb') as f:
                result_dict = pickle.load(f)
            f.close()
    

            edict_processed = {}
            for gname, gdict in result_dict['error_dict'].items():
                edict_processed[gname] = {}
                for k,v in gdict.items():
                    edict_processed[gname][(k[0],str(k[1]))] = v
            
            pm_edict = {'Gh':ExpErrorgenOp(pygsti.modelmembers.operations.LindbladErrorgen.from_elementary_errorgens(edict_processed['Gh'], state_space=QubitSpace(1))).to_dense(),
                        'Gcnot':ExpErrorgenOp(pygsti.modelmembers.operations.LindbladErrorgen.from_elementary_errorgens(edict_processed['Gcnot'], state_space=QubitSpace(2))).to_dense(),
                        'Gi':np.eye(4)}
            
            pt_dem = pt.build_pauli_twirled_model(pm_edict, circuit_to_add_noise, init_circuit=init_perfect_circuit)

            pt_distribution = compute_outcome_distribution_from_dem(pt_dem)

            with open(f'data/sweep_rounds_{rounds}_infidelity_{erate}_{i}_pauli_twirled.txt', 'wb+') as f:
                pickle.dump([pt_dem,pt_distribution], f)
            f.close()
        except:
            print(f'failed on {erate}, {i}')

    
