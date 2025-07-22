import stim
import numpy as np
import itertools
from scipy.optimize import curve_fit

format_error = lambda qubit_labels, pauli: ''.join([f'{p}{q} ' if p != 'I' else '' for p,q in zip(pauli, qubit_labels)])

decay_form = lambda x, A, b: A*(b**x)

def generate_local_depol_model(qubits, avg_strength, std):
    depol_strengths = [np.abs(np.random.normal(avg_strength, std)) for q in qubits]
    depol_layers = ''.join([f'\t DEPOLARIZE1({s}) {q}\n' for s, q in zip(depol_strengths, qubits)])
    f = np.prod([1-d for d in depol_strengths])
    return depol_strengths, depol_layers, f

def create_noisy_depolarizing_cycle(base_circuit, probs_post, errors_post, probs_pre, errors_pre, depol_layers):
    
    if len(probs_pre)==0:
        pre_meas_error = ''
    else:
        pre_meas_error = (f'\t CORRELATED_ERROR({probs_pre[0]}) {errors_pre[0]}\n')
        if len(probs_pre) > 1:
            pre_meas_error = pre_meas_error+''.join([f'\t ELSE_CORRELATED_ERROR({pr}) {err}\n' for pr, err in zip(probs_pre[1:],errors_pre[1:])])

    if len(probs_pre)==0:
        post_meas_error = ''
    else:        
        post_meas_error = (f'\t CORRELATED_ERROR({probs_post[0]}) {errors_post[0]}\n')

        if len(probs_post) > 1:
            post_meas_error = post_meas_error+''.join([f'\t ELSE_CORRELATED_ERROR({pr}) {err}\n' for pr, err in zip(probs_post[1:],errors_post[1:])])

    full_noisy_layer = pre_meas_error+depol_layers+base_circuit+post_meas_error

    return full_noisy_layer

def generate_error_model(qubits, nmeas, nrest, num_errors, num_errors_during, noise_strengths):
    #process fidelity is (1-pre_str-post_str)*(1-during_str)
    noise_strength_pre, noise_strength_during, noise_strength_post = noise_strengths
    
    valid_errors_rest = [''.join(k) for k in list(itertools.product(['I','X','Y','Z'], repeat=nrest))[1:]]
    valid_errors_meas = [''.join(k) for k in list(itertools.product(['I','X'], repeat=nmeas))[1:]]
    valid_errors_nomeas = ['I'*nmeas]
    valid_errors_pre = [''.join(k) for k in itertools.product(valid_errors_meas, valid_errors_rest)]
    valid_errors_during = [''.join(k) for k in itertools.product(valid_errors_nomeas, valid_errors_rest)]
    
    chosen_errors = np.random.choice(valid_errors_during, num_errors_during, replace=False)

    chosen_errors_pre = np.random.choice(valid_errors_pre, num_errors, replace=False)
    chosen_errors_post = np.random.choice(valid_errors_pre, num_errors, replace=False)
    #make pre and post MCM independent error
    
    errors = np.array([np.random.rand() for _ in range(num_errors)])
    errors *= noise_strength_pre/np.sum(errors)
    probs_pre = [errors[i]/(1-sum(errors[:i])) for i in range(len(errors))]
    errors_pre = [format_error(qubits, p) for p in chosen_errors_pre]
    
    errors = np.array([np.random.rand() for _ in range(num_errors)])
    errors *= noise_strength_post/np.sum(errors)
    probs_post = [errors[i]/(1-sum(errors[:i])) for i in range(len(errors))]
    errors_post = [format_error(qubits, p) for p in chosen_errors_post]

    errors = np.array([np.random.rand() for _ in range(num_errors_during)])
    errors *= noise_strength_during/np.sum(errors)
    probs_during = [errors[i]/(1-sum(errors[:i])) for i in range(len(errors))]
    errors_during = [format_error(qubits, p) for p in chosen_errors]
    
    return probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during

def generate_spam_error(qubits, prep_str, meas_str):
    prep_errors = np.array([np.random.rand() for _ in qubits])
    prep_errors *= prep_str/np.sum(prep_errors)   

    meas_errors = np.array([np.random.rand() for _ in qubits])
    meas_errors *= meas_str/np.sum(meas_errors)
    
    prep_err_layer = ''.join(f'X_ERROR({p}) {q}\n' for p,q in zip(prep_errors, qubits))
    meas_err_layer = ''.join(f'X_ERROR({p}) {q}\n' for p,q in zip(meas_errors, qubits))
    
    return prep_err_layer, meas_err_layer

def create_noisy_cycle(base_circuit, probs_post, errors_post, probs_pre, errors_pre, probs_during, errors_during):
    if len(probs_pre)==0:
        pre_meas_error = ''
    else:
        pre_meas_error = (f'\t CORRELATED_ERROR({probs_pre[0]}) {errors_pre[0]}\n')
        if len(probs_pre) > 1:
            pre_meas_error = pre_meas_error+''.join([f'\t ELSE_CORRELATED_ERROR({pr}) {err}\n' for pr, err in zip(probs_pre[1:],errors_pre[1:])])

    during_meas_error = (f'\t CORRELATED_ERROR({probs_during[0]}) {errors_during[0]}\n')
    if len(probs_pre) > 1:
        during_meas_error = during_meas_error+''.join([f'\t ELSE_CORRELATED_ERROR({pr}) {err}\n' for pr, err in zip(probs_during[1:],errors_during[1:])])
    
    if len(probs_post)==0:
        post_meas_error = ''
    else:
        post_meas_error = (f'\t CORRELATED_ERROR({probs_post[0]}) {errors_post[0]}\n')

        if len(probs_post) > 1:
            post_meas_error = post_meas_error+''.join([f'\t ELSE_CORRELATED_ERROR({pr}) {err}\n' for pr, err in zip(probs_post[1:],errors_post[1:])])

    full_noisy_layer = pre_meas_error+during_meas_error+base_circuit+post_meas_error

    return full_noisy_layer

def _select_neg_evecs(pauli, sign):
    # Selects the entries in an n-qubit that will be turned be given a -1 1Q eigenstates
    #     - pauli: The n-qubit Pauli
    #     - sign: Whether you want a -1 or +1 eigenvector
    # Returns: A bitstring whose 0/1 entries specify if you have a +1 or -1 1Q eigenstate
    n = len(pauli)
    identity_bitstring = [0 if i == 'I' else 1 for i in pauli]
    nonzero_indices = np.nonzero(identity_bitstring)[0]
    num_nid = len(nonzero_indices)
    if num_nid % 2 == 0:
        if sign == 1:
            choices = np.arange(start = 0, stop = num_nid+1, step = 2)
        else:
            choices = np.arange(start = 1, stop = num_nid, step = 2)
    else:
        if sign == 1:
            choices = np.arange(start = 0, stop = num_nid, step = 2)
        else:
            choices = np.arange(start = 1, stop = num_nid+1, step = 2)
    num_neg_evecs = np.random.choice(choices)
    assert((-1)**num_neg_evecs == sign)
    neg_evecs = np.random.choice(nonzero_indices, num_neg_evecs, replace = False)
    
    bit_evecs = np.zeros(n)
    bit_evecs[neg_evecs] = 1
    assert('I' not in np.array(pauli)[nonzero_indices])
    
    return bit_evecs

def generate_circuits(prep_pauli, depth, n_circs, qubits, full_noisy_layer, prep_err='', meas_err=''):
    clist = []
    signs = []
    meas= ''.join(f'M{pauli} {q}\n' if pauli !='I' else f'MZ {q}\n' for pauli, q in zip(prep_pauli, qubits))
    core = f'REPEAT {depth}{{\n'+full_noisy_layer+'}\n'    
    
    for _ in range(n_circs):
        if not all([p=='I' for p in prep_pauli]):
            sign = np.random.choice([-1,1])
        else:
            sign = 1
        prep = make_prep_layers(prep_pauli, sign, qubits)
        if depth>0:
            #print(prep+prep_err+core+meas_err+meas)
            clist.append(stim.Circuit(prep+prep_err+core+meas_err+meas))
        else:
            
            clist.append(stim.Circuit(prep+prep_err+meas_err+meas))
        signs.append(sign)
    
    return clist, signs

def make_prep_layers(prep_pauli, sign, qubits):
    basis_gate = {'X':'H {}\n','Y':'H_YZ {}\n'}
    prep_layer = ''.join([basis_gate[p].format(q) if p=='X' or p=='Y' else '' for p,q in zip(prep_pauli, qubits)])
    #TODO: add sign-changing layer
    sign_bits = _select_neg_evecs([p for p in prep_pauli], sign)
    sign_layer = ''.join([f'X {q}\n' if b else '' for b,q in zip(sign_bits,qubits)])

    return sign_layer+prep_layer

def compute_avg_val(datalist, prep_pauli, toggle, d, signs):
    energies = []
    final_meas = [True if p != 'I' else False for p in prep_pauli]
    target_string = toggle*d+final_meas
    #need to filter out anything that isn't toggled within MCMs
    #for mult qubit measurement, toggle is a list of booleans
    for counts, sign in zip(datalist,signs):
        products = [sign*(-1)**int(np.sum(np.logical_and(r, target_string))) for r in counts]
        energy = sum(products)/len(products)
        energies.append(energy)
    
    return energies, np.mean(energies)

def estimate_process_fidelity(k, nrest, num_meas, cycle, num_randomizations, depths, prep_err, meas_err, nshots=1000, all_paulis = False, show_progress=True):
    paulis=list(itertools.product(['I','X','Y','Z'], repeat=nrest))
    
    toggle_paulis = list(itertools.product([False, True], repeat=num_meas))
    qubits = range(nrest+num_meas)
    n = nrest+num_meas
    if all_paulis:
        assert(n < 7)
        #generate all Paulis
        random_paulis=list(itertools.product(['I','X','Y','Z'], repeat=nrest))
    else: 
        #should probably generate on the fly
        inds = np.random.choice(len(paulis), size=k, replace=False)
        random_paulis = [paulis[ind] for ind in inds]
    
    #full Paulis: I and Z on measured qubits
    measured_paulis = list(itertools.product(['I','Z'],repeat=num_meas))
    required_paulis = list(sum(j,()) for j in itertools.product(measured_paulis, random_paulis))
    
    estimates = []
    estimates_std = []
    evals_by_p = {}
    tvals_by_p = {}
    all_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
    avg_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
    for pauli_tuple in required_paulis:
        pauli = ''.join(pauli_tuple)
        

        if show_progress:
            print(pauli)
        #generate and run experiments for these Paulis
        for depth in depths:
            circuits, signs = generate_circuits(pauli, depth, num_randomizations, qubits, cycle, prep_err, meas_err)

            datalist = [c.compile_sampler().sample(shots=nshots) for c in circuits]
            

            for toggle in toggle_paulis:
                energies, avg_energy = compute_avg_val(datalist, pauli, list(toggle), depth, signs)
                avg_energies[(pauli, toggle)].append(avg_energy)
                all_energies[(pauli,toggle)].append(energies)

        for toggle in toggle_paulis:
            res = curve_fit(decay_form, depths, avg_energies[(pauli,toggle)], p0=[1,0.99])[0]
            lambda_est = res[-1]
            estimates.append(lambda_est)
            
            #if we're also sampling Paulis, these should be resampled in the bootstrap.
            est_std = bootstrap_error(all_energies[(pauli,toggle)], depths, bootstrap_samples=100, n_shots=nshots)
            estimates_std.append(est_std)

            evals_by_p[(pauli, toggle)] = lambda_est

    est = sum(estimates)/len(estimates)

    est_std = np.mean(estimates_std)
    return est, est_std, evals_by_p, estimates_std, all_energies

def bootstrap_error(data_by_depth, depths, bootstrap_samples=100, n_shots=1000):
    #actually what should I resample here? 
    #if I were actually RCing, I could resample by circuit, but here there's only one circuit to run. 
    rs = []
    for _ in range(bootstrap_samples):
        new_ps = {}
        for energies, d in zip(data_by_depth,depths):
            if type(energies) != list:
                energies = [energies]
            probs = [(e+1)/2 for e in energies]
            #sample probabilities from existing
            sampled_probs = np.random.choice(probs, len(energies))

            #from this set of probs, take samples from a binomial distribution
            new_ps[d] = []
            for p in sampled_probs:
                outcomes = np.random.binomial(1, p, size=n_shots)
                #compute new  probs
                new_p = sum(outcomes)/len(outcomes)
                new_ps[d].append(new_p)
        
        mean_sps = [2*np.mean(new_ps[d])-1 for d in depths]
        result = curve_fit(decay_form, [d for d in depths if d != 1], mean_sps, p0=[1,0.99])[0]
        rs.append(result[-1])
    return np.std(rs)


def bootstrap_error_fidelity(data_by_depth_by_pauli, depths, paulis, toggles, bootstrap_samples=100, n_shots=1000): 
    rs = []
    for _ in range(bootstrap_samples):
        sampled_paulis = np.random.choice(paulis, len(paulis))
        for pauli in sampled_paulis:
            for toggle in toggles:
                data_by_depth = data_by_depth_by_pauli[(pauli, toggle)]
                new_ps = {}
                for energies, d in zip(data_by_depth,depths):
                    if type(energies) != list:
                        energies = [energies]
                    probs = [(e+1)/2 for e in energies]
                    #sample probabilities from existing
                    sampled_probs = np.random.choice(probs, len(energies))
                    new_ps[d] = []
                    for p in sampled_probs:
                        outcomes = np.random.binomial(1, p, size=n_shots)
                        #compute new  probs
                        new_p = sum(outcomes)/len(outcomes)
                        new_ps[d].append(new_p)
        
        mean_sps = [2*np.mean(new_ps[d])-1 for d in depths]
        result = curve_fit(decay_form, [d for d in depths if d != 1], mean_sps, p0=[1,0.99])[0]
        rs.append(result[-1])
    return np.std(rs)


def estimate_process_fidelity_alt(k, nrest, num_meas, cycle, num_randomizations, depths, prep_err, meas_err, nshots=1000, all_paulis = False, show_progress=True):
    unmeas_paulis=list(itertools.product(['I','X','Y','Z'], repeat=nrest))
    measured_paulis = list(itertools.product(['I','Z'],repeat=num_meas))
    paulis = [p2+p1 for p1, p2 in list(itertools.product(unmeas_paulis, measured_paulis))]
    
    toggle_paulis = list(itertools.product([False, True], repeat=num_meas))
    qubits = range(nrest+num_meas)
    if all_paulis:
        assert(n < 7)
        random_paulis=list(itertools.product(['I','X','Y','Z'], repeat=nrest))
    else: 

        inds = np.random.choice(len(paulis), size=k, replace=False)
        random_paulis = [paulis[ind] for ind in inds]
        inds = np.random.choice(len(toggle_paulis), size=k, replace=True)
        random_toggles = [toggle_paulis[ind] for ind in inds]

    required_paulis = random_paulis
    estimates = []
    estimates_std = []
    evals_by_p = {}
    tvals_by_p = {}
    all_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
    avg_energies = {(''.join(p for p in pauli),toggle):[] for toggle in toggle_paulis for pauli in required_paulis}
    for pauli_tuple, toggle in zip(required_paulis,random_toggles):
        pauli = ''.join(pauli_tuple)
        

        if show_progress:
            print(pauli)
        #generate and run experiments for these Paulis
        for depth in depths:
            circuits, signs = generate_circuits(pauli, depth, num_randomizations, qubits, cycle, prep_err, meas_err)

            datalist = [c.compile_sampler().sample(shots=nshots) for c in circuits]

            energies, avg_energy = compute_avg_val(datalist, pauli, list(toggle), depth, signs)
            avg_energies[(pauli, toggle)].append(avg_energy)
            all_energies[(pauli,toggle)].append(energies)

        res = curve_fit(decay_form, depths, avg_energies[(pauli,toggle)], p0=[1,0.99])[0]
        lambda_est = res[-1]
        estimates.append(lambda_est)
            
            #if we're also sampling Paulis, these should be resampled in the bootstrap.
        est_std = bootstrap_error(all_energies[(pauli,toggle)], depths, bootstrap_samples=100, n_shots=nshots)
        estimates_std.append(est_std)

        evals_by_p[(pauli, toggle)] = lambda_est #, est_std)
            
        
    #Take sum
    est = sum(estimates)/len(estimates)
    est_std = np.mean(estimates_std)
    return est, est_std, evals_by_p, estimates_std, all_energies