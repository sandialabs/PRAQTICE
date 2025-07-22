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

import stim
import numpy as np
from mcmcb_stim_tools import * 

import pickle

param_list= np.arange(0.0001, 0.0201, 0.00025)
#pre, during, post
noise_strengths = [0.02, 0.04, 0.02]

rs = [0.1, 1, 10]

num_randomizations = 30
depths = [0,2,4,8,16,32]


unmeas_widths = [4,6]
meas_widths = [1,2]

if __name__=="__main__":
    np.random.seed(rank)
 
    for offset in [0,1]:
        idx = rank+40*offset
        if rank==0:
            print(idx)

        for r in range(1):
            for cpd in range(1):
                for nrest in unmeas_widths:
                    for nmeas in meas_widths:
                        if rank==0:
                            print(nmeas, nrest)
                        
                        k = np.min([4**nrest,100])
                        if rank==0:
                            print(k)
                        
                        n=nmeas+nrest
                        
                        qubits = range(n)
                        
                        
                        prep_str = param_list[idx]*n/2
                        
                        meas_str = param_list[idx]*n/2
                        
                        if rank==0:
                            print(prep_str, meas_str)
                        
                        base_circuit = ''.join(f'\t M {q}\n' for q in qubits[:nmeas])
                        
                        noise_strengths = [0.01, 0.02, 0.01]
                        
                        num_errors = 4**nrest-1
                        
                        probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during = generate_error_model(qubits, nmeas, nrest, num_errors, num_errors, noise_strengths)
                        
                        full_noisy_layer = create_noisy_cycle(base_circuit, probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during)

                        prep_err, meas_err = generate_spam_error(qubits, prep_str, meas_str)


                        if rank==0:
                            prog = True
                        else:
                            prog = False
                        est, est_std, evals_by_p, stds, energies = estimate_process_fidelity_alt(k, nrest, nmeas, full_noisy_layer,num_randomizations, depths, prep_err, meas_err, show_progress=prog, nshots=10000)                   


                        result_dict = {'probs_post':probs_post, 'errors_post':errors_post, 'probs_pre':probs_pre, 
                                       'errors_pre':errors_pre, 'probs_during':probs_during, 'errors_during':errors_during, 
                                       'fidelity':0.99*0.98**2, 'fidelity_est':est, 'fidelity_std':est_std,
                                       'evals_by_p':evals_by_p, 'energies':energies, 'stds':stds}

                        with open(f'results/stim/sweep_spam/{nmeas}_{nrest}_qubit_mcmcb_run_{idx}_semi_dense.pkl', 'wb+') as f:
                            pickle.dump(result_dict, f)
