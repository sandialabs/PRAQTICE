import numpy as np
import pygsti

def sample_error_rates(strengths, n):
    '''
    Samples an error rates dictionary for dependent gates.
    '''

    error_rates_dict = {}
    # Sample stochastic error rates. First we sample the overall stochastic error rate.
    # Then we sample (and normalize) the individual stochastic error rates
    stochastic_strength = strengths['S'] * np.random.random()
    s_error_rates = np.random.random(4 ** n - 1)
    s_error_rates = s_error_rates / np.sum(s_error_rates) * stochastic_strength

    hamiltonian_strength = strengths['H'] * np.random.random()
    h_error_rates = np.random.random(4 ** n - 1)
    h_error_rates = h_error_rates * np.sqrt(hamiltonian_strength) / np.sqrt(np.sum(h_error_rates**2))

    error_rates_dict.update({('S', i + 1): s_error_rates[i] for i in range(4 ** n - 1)})
    error_rates_dict.update({('H', i + 1): h_error_rates[i] for i in range(4 ** n - 1)})

    return error_rates_dict

def sample_error_rates_dict(pspec, strengths):
    """
    For example:
        strengths = {1: {'S':0.001, 'H':0.01}, 
                     2: {'S':0.01,'H':0.1}}

    The 'S' and 'H' entries in the strengths dictionary give the maximum possible contribution to the infidelity from a given gate.
    """
    qubits = pspec.qubit_labels
    errors_rates_dict = {}
    for gate, availability in pspec.availability.items():
        n = pspec.gate_num_qubits(gate)
        if availability == 'all-edges':
            assert(n == 1), "Currently require all 2-qubit gates have a specified availability!"
            qubits_for_gate = qubits
        else:
            qubits_for_gate = availability  
        for qs in qubits_for_gate:
            label = pygsti.baseobjs.Label(gate, qs)
            # First, check if there's a strength specified for this specific gate.
            max_stength = strengths.get(label, None) # to get highly biased errors can set generic error rates to be low, then set it to be high for one or two particular gates.
            # Next, check if there's a strength specified for all gates with this name
            if max_stength is None:
                max_stength = strengths.get(gate, None)
            # Finally, get error rate for all gates on this number of qubits.
            if max_stength is None:
                max_stength = strengths[n]
            # Sample error rates.
            errors_rates_dict[label] = sample_error_rates(max_stength, n)
 
    return errors_rates_dict

def construct_model(pspec, errors_dict):
    lindblad_error_coeffs = {}
    for label in errors_dict:
        lindblad_error_coeffs[label] = errors_dict[label]
    return pygsti.models.create_crosstalk_free_model(pspec, lindblad_error_coeffs=lindblad_error_coeffs, 
                                                     custom_gates=custom_gates)
