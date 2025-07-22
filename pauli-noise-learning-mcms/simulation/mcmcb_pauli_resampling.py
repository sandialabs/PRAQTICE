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


inf = 0.02
avg_strength=0.01
std = 0.002

params = [[0.01, 0.04, 0.01], [0, 0.01, 0.05],[0.03, 0, 0.03]]

nrest = 5 #1024 Paulis
num_meas = 2 #4 Paulis

#plan to do 32 processes, each running 32 experiments
nshots=1000

#k = 20
num_randomizations = 20
depths = [0,2,4,8,16,32]

if __name__=="__main__":
    np.random.seed(rank)

    for mnum, plist in enumerate(params):

        noise_strengths = plist
        n=num_meas+nrest

        qubits = range(n)

        num_errors = 2**n

        prep_str = 0.005*n

        meas_str = 0.01*n

        base_circuit = ''.join(f'\t M {q}\n' for q in qubits[:num_meas])

        #depol_strengths, depol_layers, f = generate_local_depol_model(qubits[num_meas:], avg_strength, std)

        probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during = generate_error_model(qubits, num_meas, nrest, num_errors, num_errors, noise_strengths)

        full_noisy_layer = create_noisy_cycle(base_circuit, probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during)

        #full_noisy_layer = create_noisy_depolarizing_cycle(base_circuit, probs_post, errors_post, probs_pre, errors_pre, depol_layers)

        prep_err, meas_err = generate_spam_error(qubits, prep_str, meas_str)


        paulis=list(itertools.product(['I','X','Y','Z'], repeat=nrest))

        toggle_paulis = list(itertools.product([False, True], repeat=num_meas))
        #print(list(toggle_paulis))
        qubits = range(nrest+num_meas)

        random_paulis = paulis[32*rank:32*(rank+1)]

        #full Paulis: I and Z on measured qubits
        measured_paulis = list(itertools.product(['I','Z'],repeat=num_meas))
        required_paulis = list(sum(j,()) for j in itertools.product(measured_paulis, random_paulis))

        estimates = []
        estimates_std = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
        evals_by_p = {}
        tvals_by_p = {}
        all_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
        avg_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}

        for j, pauli_tuple in enumerate(required_paulis):
            pauli = ''.join(pauli_tuple)
            if rank==0: print(j)

            #generate and run experiments for these Paulis
            for depth in depths:
                #print(depth)
                circuits, signs = generate_circuits(pauli, depth, num_randomizations, qubits, full_noisy_layer, prep_err, meas_err)

                datalist = [c.compile_sampler().sample(shots=nshots) for c in circuits]


                for toggle in toggle_paulis:
                    energies, avg_energy = compute_avg_val(datalist, pauli, list(toggle), depth, signs)
                    avg_energies[(pauli, toggle)].append(avg_energy)
                    all_energies[(pauli,toggle)].append(energies)

            for toggle in toggle_paulis:
                res = curve_fit(decay_form, depths, avg_energies[(pauli,toggle)], p0=[1,0.99])[0]
                lambda_est = res[-1]
                estimates.append(lambda_est)

                est_std = bootstrap_error(all_energies[(pauli,toggle)], depths, bootstrap_samples=100, n_shots=nshots)
                estimates_std[(pauli, toggle)] = est_std

                evals_by_p[(pauli, toggle)] = lambda_est #, est_std)

        est = sum(estimates)/len(estimates)             


        result_dict = {'probs_post':probs_post, 'errors_post':errors_post, 'probs_pre':probs_pre, 
                       'errors_pre':errors_pre, 'probs_during':probs_during, 'errors_during':errors_during, 
                       'fidelity':(1-plist[1])*(1-plist[0])*(1-plist[2]), 'fidelity_est':est,
                       'evals_by_p':evals_by_p, 'estimate_stds': estimates_std}

        with open(f'results/all-paulis/{num_meas}_{nrest}_qubit_mcmcb_run_{rank}_plist_{mnum}.pkl', 'wb+') as f:
            pickle.dump(result_dict, f)
        
    inf = 0.02
    avg_strength=0.01
    std = 0.002

    params = [[0.01, 0.01, 0.01],[0, 0.005,0.02]]

    for mnum, plist in enumerate(params):

        noise_strengths = plist
        n=num_meas+nrest

        qubits = range(n)

        num_errors = 2**n

        prep_str = 0.005*n

        meas_str = 0.01*n

        base_circuit = ''.join(f'\t M {q}\n' for q in qubits[:num_meas])

        depol_strengths, depol_layers, f = generate_local_depol_model(qubits[num_meas:], avg_strength, std)

        probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during = generate_error_model(qubits, num_meas, nrest, num_errors, num_errors, noise_strengths)


        full_noisy_layer = create_noisy_depolarizing_cycle(base_circuit, probs_post, errors_post, probs_pre, errors_pre, depol_layers)

        prep_err, meas_err = generate_spam_error(qubits, prep_str, meas_str)


        paulis=list(itertools.product(['I','X','Y','Z'], repeat=nrest))

        toggle_paulis = list(itertools.product([False, True], repeat=num_meas))
        qubits = range(nrest+num_meas)

        random_paulis = paulis[32*rank:32*(rank+1)]

        #full Paulis: I and Z on measured qubits
        measured_paulis = list(itertools.product(['I','Z'],repeat=num_meas))
        required_paulis = list(sum(j,()) for j in itertools.product(measured_paulis, random_paulis))

        estimates = []
        estimates_std = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
        evals_by_p = {}
        tvals_by_p = {}
        all_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
        avg_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}

        for j, pauli_tuple in enumerate(required_paulis):
            pauli = ''.join(pauli_tuple)
            if rank==0: print(j)

            #generate and run experiments for these Paulis
            for depth in depths:
                circuits, signs = generate_circuits(pauli, depth, num_randomizations, qubits, full_noisy_layer, prep_err, meas_err)

                datalist = [c.compile_sampler().sample(shots=nshots) for c in circuits]

                for toggle in toggle_paulis:
                    energies, avg_energy = compute_avg_val(datalist, pauli, list(toggle), depth, signs)
                    avg_energies[(pauli, toggle)].append(avg_energy)
                    all_energies[(pauli,toggle)].append(energies)

            for toggle in toggle_paulis:
                res = curve_fit(decay_form, depths, avg_energies[(pauli,toggle)], p0=[1,0.99])[0]
                lambda_est = res[-1]
                estimates.append(lambda_est)

                est_std = bootstrap_error(all_energies[(pauli,toggle)], depths, bootstrap_samples=100, n_shots=nshots)
                estimates_std[(pauli, toggle)] = est_std

                evals_by_p[(pauli, toggle)] = lambda_est

        est = sum(estimates)/len(estimates)

        result_dict = {'probs_post':probs_post, 'errors_post':errors_post, 'probs_pre':probs_pre, 
                       'errors_pre':errors_pre, 'probs_during':probs_during, 'errors_during':errors_during, 
                       'fidelity':f*(1-plist[0])*(1-plist[2]), 'fidelity_est':est,
                       'evals_by_p':evals_by_p, 'estimate_stds': estimates_std}

        with open(f'results/all-paulis/{num_meas}_{nrest}_qubit_mcmcb_run_{rank}_depol_{mnum}.pkl', 'wb+') as f:
            pickle.dump(result_dict, f)
