import numpy as np
import pygsti
import pygsti.baseobjs.label as _lbl
import pygsti.modelmembers as mm
from pygsti.circuits import Circuit as _cir
import pygsti.algorithms.randomcircuit as _rc
from pygsti.tools import symplectic as _symp
import itertools
from pygsti.processors import CliffordCompilationRules as CCR
from scipy.optimize import curve_fit

def get_clifford_from_unitary(U):
    clifford_unitaries = {k: v for k, v in pygsti.tools.internalgates.standard_gatename_unitaries().items()
                      if 'Gc' in k and v.shape == (2, 2)}
    
    for k,v in clifford_unitaries.items():
        for phase in [1, -1, 1j, -1j]:
            if np.allclose(U, phase*v):
                return k
    
    raise RuntimeError(f'Failed to look up Clifford for unitary:\n{U}')

    
def pauli_vector_to_clifford_layer(p, qubits):
    
    n = len(qubits)
    layer = []
    for i, q in enumerate(qubits):

        if p[i] == 0 and p[i+n] == 0:  # I
            label = 'Gc0'
        elif p[i] == 2 and p[i+n] == 0:  # Z
            label = 'Gc9'
        elif p[i] == 0 and p[i+n] == 2:  # X
            label = 'Gc3'
        elif p[i] == 2 and p[i+n] == 2:  # Y
            label = 'Gc6'

        layer.append(pygsti.baseobjs.Label(label, q))

    return pygsti.baseobjs.Label(layer)

clifford_label_cache = {}
def update_clifford_labels(layer, p, q, qubits):
    used_qubits = []

    new_layer = []
    n = len(qubits)

    clifford_unitaries = {k: v for k, v in pygsti.tools.internalgates.standard_gatename_unitaries().items()
                      if 'Gc' in k and v.shape == (2, 2)}

    Ps = pauli_vector_to_clifford_layer(p, qubits)
    Qs = pauli_vector_to_clifford_layer(q, qubits)

    for g in layer:
        idx = qubits.index(g.qubits[0])
        p = Ps[idx]
        q = Qs[idx]

        cPUQ = clifford_label_cache.get((p.name, g.name, q.name), None)

        if cPUQ is None:
            # Get unitary
            U = clifford_unitaries[g.name]
            P = clifford_unitaries[p.name]
            Q = clifford_unitaries[q.name]

            # Do RC
            PUQ = Q @ U @ P

            # Look up RC'd clifford (up to phase)
            cPUQ = get_clifford_from_unitary(PUQ)


            # Store in cache for efficiency
            clifford_label_cache[p.name, g.name, q.name] = cPUQ

        new_label = pygsti.baseobjs.Label(cPUQ, g.qubits[0])
        new_layer.append(new_label)
        used_qubits.append(g.qubits[0])

    new_layer = pygsti.baseobjs.Label(new_layer)
    assert(set(used_qubits) == set(qubits))

    return new_layer

def pad_layer_clifford(layer, qubits):

    padded_layer = list(layer)
    used_qubits = []
    for g in layer:
        for q in g.qubits:
            used_qubits.append(q)

    for q in qubits:
        if q not in used_qubits:
            padded_layer.append(pygsti.baseobjs.Label('Gc0', q))

    return pygsti.baseobjs.Label(padded_layer)

def pauli_randomize_clifford_circuit(circ):

    d = circ.depth
    n = circ.width
    p = np.zeros(2*n, int)
    q = np.zeros(2*n, int)
    rc_circ = pygsti.circuits.Circuit(line_labels=circ.line_labels, editable=True)
    qubits = circ.line_labels
    
    measurement_bitflips = []

    for i in range(d):

        layer = circ.layer_label(i).components
        

        if len(layer) == 0 or layer[0].name in [f'Gc{i}' for i in range(24)]:
            q = 2 * np.random.randint(0, 2, 2*n)
            padded_layer = pad_layer_clifford(layer, qubits)
            rc_layer = update_clifford_labels(padded_layer, p, q, qubits)
            rc_circ.insert_layer_inplace(rc_layer, i)
            p = q
            
        else:
            rc_circ.insert_layer_inplace(layer, i)
            for g in layer:
                if g.name == 'Gcnot':
                    (control, target) = g.qubits
                    p[qubits.index(control)] = (p[qubits.index(control)] + p[qubits.index(target)]) % 4
                    p[n + qubits.index(target)] = (p[n + qubits.index(control)] + p[n + qubits.index(target)]) % 4
                if g.name == 'Gcphase':
                    (control, target) = g.qubits
                    #target - forward propagate phase based on state of control
                    p[qubits.index(target)] = (p[n+qubits.index(control)] + p[qubits.index(target)]) % 4
                    #control- back propagate a phase based on state of target
                    p[qubits.index(control)] = (p[qubits.index(control)] + p[n + qubits.index(target)]) % 4
                    
                elif g.name=='Gc0':
                    pass
                elif g.name=='Iz':    
                    idx = qubits.index(g.qubits[0])
                    #need to record if X was applied pre-measurement
                    x = int(p[idx+n]==2)
                    measurement_bitflips.append(str(x))
                else:
                    raise ValueError("Circuit can only contain Gcnot, Gcphase, and Gc[0-23] gates in separate layers!")

    bs = ''.join(measurement_bitflips+[str(b // 2) for b in q[n:]])

    rc_circ.done_editing()

    return rc_circ, bs, q

def generate_circuits_rc(cycle, depths, required_paulis, compilations, num_randomizations, qubit_labels, show_progress=True):  
    cs_by_pauli = {}
    signs_by_pauli= {}
    tbs_by_pauli = {}
    for p in required_paulis:
        if show_progress:
            print(p)
        clist = []
        signlist = []
        tbslist = []
        for d in depths:
            cs = []
            signs = []
            tbss = []
            for _ in range(num_randomizations):
                if not all([i=='I' for i in p]):
                    sign = np.random.choice([-1,1])
                else: 
                    sign = 1
                stab_state, stab_phase, s_prep, p_prep, prep_layer = _rc._sample_stabilizer(p, sign, compilations['absolute'], qubit_labels)
                s_layer, p_layer, meas_layer = _rc._stabilizer_to_all_zs(p, qubit_labels, compilations['absolute'])
                circuit = prep_layer.serialize().copy(editable=True)
                circuit.append_circuit_inplace(cycle.repeat(d))
                circuit.append_circuit_inplace(meas_layer.serialize())
                circuit.delete_idle_layers_inplace()
                circuit.done_editing()
                circuit, tbs, _ = pauli_randomize_clifford_circuit(circuit)
                cs.append(circuit)
                tbss.append(tbs)
                signs.append(sign)
            clist.append(cs)
            signlist.append(signs)
            tbslist.append(tbss)
                #circuit is prep-(cycle)^d-meas
        cs_by_pauli[p] = clist
        signs_by_pauli[p] = signlist
        tbs_by_pauli[p] = tbslist

    return cs_by_pauli, signs_by_pauli, tbs_by_pauli

###Analysis###

def determine_new_sign(c, tbs, pspec, meas_pauli):
    n = c.width
    #remove measured qubit
    rest_c = c.copy(editable=True)
    rest_c = rest_c.replace_gatename_with_idle('Iz')
    init_stab, init_phase = _symp.prep_stabilizer_state(len(rest_c.line_labels))
    s_c, p_c = _symp.symplectic_rep_of_clifford_circuit(rest_c, pspec=pspec.subset(gate_names_to_include='all', qubit_labels_to_keep=rest_c.line_labels)) 
    s_state, p_state = _symp.compose_cliffords(s1 = init_stab, p1 = init_phase, s2 = s_c, p2 = p_c)
    #measure at end
    meas = _rc._measure(s_state, p_state) #list of 0 and 1 outcomes
    meas_to_count = [0 if p=='I' else 1 for p in meas_pauli]
    new_sign = (-1)**(np.dot(meas, meas_to_count))
    new_tbs = tbs[:-n]+'0'*n
    return new_sign, new_tbs

def outcome_energy(outcome, measurement, sign, tbs):
    energy = 1
    for i,j,k in zip(outcome,measurement,tbs):
        if j == 'Z' and (i != k):
            energy = -1*energy
    return sign*energy

def avg_energy(cd, measurement, sign, tbs):
    energy = 0
    total = sum(cd.values())
    for i,count in cd.items():
        out_eng = outcome_energy(i[0],measurement,sign, tbs)
        energy += count * out_eng    
    return energy / total

def ignore_mcm_results(count_dict):
    new_dict = {}
    for key, count in count_dict.items():
        new_key = (key[-1],)
        if new_key in new_dict:
            new_dict[new_key] += count
        else:
            new_dict[new_key] = count
    return new_dict

decay_form = lambda x, A, b: A*(b**x)

#needs to be generalized to multiple qubits
def avg_energy_sign_mod(cd, measurement, sign, tbs, measured_qs, toggled_qs):
    energy = 0
    total = sum(cd.values())
    for i,count in cd.items():
        if len(i)>1:
            #is this right tbs
            if len(toggled_qs)>0:
                mcm_results = [0 if b=='p0' else 1 for b in i[:-1]]
                counted_mcm_results = [1 if q in toggled_qs else 0 for q in measured_qs]*(len(mcm_results)//len(measured_qs))
                mcm_tbs = [0 if tbs[j]=='0' else 1 for j in range(len(mcm_results))]
                mcm_adjusted_results = np.dot(np.logical_xor(mcm_results, mcm_tbs), counted_mcm_results)
                this_sign = (-1)**(mcm_adjusted_results%2)
            
            else:
                this_sign=1
        else:
            this_sign=1

        out_eng = outcome_energy(i[-1],measurement,sign, tbs[-len(measurement):])*this_sign
        energy += count * out_eng    
    return energy / total

def compute_eigenvalue_decays(data_by_pauli, cs_by_pauli, signs_by_pauli, tbs_by_pauli):
    energies_by_pauli = {}
    circuit_energies_by_pauli = {}
    for pauli, ds_by_d in data_by_pauli.items():
        circuits = cs_by_pauli[pauli]
        signs = signs_by_pauli[pauli]
        tbs = tbs_by_pauli[pauli]
        energies = []
        #transform into z type Pauli
        meas_pauli = [p if p in ['I', 'Z'] else 'Z' for p in pauli]
        circuit_energies_by_pauli[pauli] = []
        avg_energies = []
        for clist, signlist, ds in zip(circuits, signs, ds_by_d):
            circuit_energies = []
            
            for c, sign in zip(clist,signlist):
                dsrow = ds[c]
                cd = ignore_mcm_results(dsrow.to_dict())
                energy = avg_energy(cd, meas_pauli, sign)
                circuit_energies.append(energy)
            avg_energies.append(np.mean(circuit_energies))
            circuit_energies_by_pauli[pauli].append(circuit_energies)
        energies_by_pauli[pauli] = avg_energies
            
    return energies_by_pauli, circuit_energies_by_pauli

def compute_toggle_decays(data_by_pauli, cs_by_pauli, signs_by_pauli, tbs_by_pauli):
    #compute pauli measurement results, use MCM results
    energies_by_pauli = {}
    circuit_energies_by_pauli = {}
    for pauli, ds_by_d in data_by_pauli.items():
        circuits = cs_by_pauli[pauli]
        signs = signs_by_pauli[pauli]
        energies = []
        #transform into z type Pauli
        circuit_energies_by_pauli[pauli] = []
        meas_pauli = [p if p in ['I', 'Z'] else 'Z' for p in pauli]
        avg_energies = []
        for clist, signlist, ds in zip(circuits, signs, ds_by_d):
            circuit_energies = []
            
            for c, sign in zip(clist,signlist):
                meas_pauli = [p if p in ['I', 'Z'] else 'Z' for p in pauli]
                dsrow = ds[c]
                cd = ignore_mcm_results(dsrow.to_dict())
                energy = avg_energy_sign_mod(dsrow.to_dict(), meas_pauli, sign)
                circuit_energies.append(energy)
            avg_energies.append(np.mean(circuit_energies))
            circuit_energies_by_pauli[pauli].append(circuit_energies)
        energies_by_pauli[pauli] = avg_energies
    return energies_by_pauli, circuit_energies_by_pauli

def bootstrap_error(data_by_depth, depths, bootstrap_samples=100, n_shots=10000): 
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
        result = curve_fit(decay_form, depths, mean_sps, p0=[1,0.99])[0]
        rs.append(result[-1])
    return np.std(rs)