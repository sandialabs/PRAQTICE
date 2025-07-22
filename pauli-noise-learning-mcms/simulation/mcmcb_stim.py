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



#parameters
#num_errors = 2
#num_errors_during= 2

param_list= np.arange(0.0001, 0.0601, 0.0005)
#pre, during, post
noise_strengths = [0.02, 0.01, 0.02]



#k = 20
num_randomizations = 10
depths = [0,2,4,8,16,32]


unmeas_widths = [1,2,4,6,8]
meas_widths = [1,2]

if __name__=="__main__":
    np.random.seed(rank)
 
    for offset in [0,1,2]:
        idx = rank+40*offset
        if rank==0:
            print(idx)

        for r in range(1):
            for cpd in range(1):
                for nrest in unmeas_widths:
                    for nmeas in meas_widths:
                        if rank==0:
                            print(nmeas, nrest)
                        
                        inf=param_list[idx]
                        k = np.min([4**nrest,100])
                        if rank==0:
                            print(k)
                        
                        noise_strengths = [0.5*inf, inf, 0.5*inf]
                        n=nmeas+nrest
                        
                        qubits = range(n)
                        
                        num_errors = 3**nrest
                        
                        prep_str = 0.005*n
                        
                        meas_str = 0.01*n
                        
                        base_circuit = ''.join(f'\t M {q}\n' for q in qubits[:nmeas])
                        
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
                                       'fidelity':1-inf, 'fidelity_est':est, 'fidelity_std':est_std,
                                       'evals_by_p':evals_by_p, 'energies':energies, 'stds':stds}

                        with open(f'results/stim/{nmeas}_{nrest}_qubit_mcmcb_run_{idx}_alt.pkl', 'wb+') as f:
                            pickle.dump(result_dict, f)
