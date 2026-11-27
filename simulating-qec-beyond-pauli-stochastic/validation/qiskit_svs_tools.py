import pygsti
import numpy as np
import stim
from qiskit.quantum_info import SparsePauliOp,StabilizerState, Statevector
from qiskit.primitives import StatevectorEstimator
from qiskit import QuantumCircuit

from qiskit.primitives import StatevectorSampler
from collections import Counter

def build_rep_code_qiskit(d):
    qiskit_circ=QuantumCircuit(2*d-1)
    #even numbers: data
    #odd numbers: ancilla
    dqs = [2*i for i in range(d)]
    aqs = [2*i+1 for i in range(d-1)]
    for aq in aqs:
        qiskit_circ.cx(aq-1,aq)
        qiskit_circ.cx(aq+1,aq)
    return qiskit_circ

# we need the different error gates to commute for this to perfectly correspond to our true model
#for now I'm just going to assume that to avoid figuring out fully general SU(4) gates in Qiskit
def create_noisy_qiskit_rep_code_circuit(d, error_rates_dict, rounds=1):
    #TODO make error model same as in approximate sims
    qiskit_circ=QuantumCircuit(2*d-1+(d-1)*(rounds-1))
    #even numbers: data
    #odd numbers: ancilla
    dqs = [2*i for i in range(d)]
    aqs = [2*i+1 for i in range(d-1)]
    for r in range(rounds):
        for aq in aqs:
            #compute virtual qubits
            if r==0:
                vaq = aq
            elif r==1:
                vaq = 2*d-1+((aq-1)//2)
            else:
                vaq = 2*d-1+(r-1)*(d-1)+((aq-1)//2)
            qiskit_circ.cx(aq-1,vaq)
            for k,v in error_rates_dict.items():
                if k[0]==0:
                    qiskit_circ = append_error_qiskit_circ(qiskit_circ, pauli, v)
            qiskit_circ.rx(2*error_rates_dict['Gcnot','IX'],vaq)           
            qiskit_circ.rx(2*error_rates_dict['Gcnot','XI'],aq-1)
            qiskit_circ.rxx(2*error_rates_dict['Gcnot','XX'],vaq,aq-1)
            
        for aq in aqs:
            #compute virtual qubits
            if r==0:
                vaq = aq
            elif r==1:
                vaq = 2*d-1+((aq-1)//2)
            else:
                vaq = 2*d-1+(r-1)*(d-1)+((aq-1)//2)
            for k,v in error_rates_dict.items():
                if k[0]==0:
                    qiskit_circ = append_error_qiskit_circ(qiskit_circ, pauli, v)
            qiskit_circ.cx(aq+1,vaq)
            qiskit_circ.rx(2*error_rates_dict['Gcnot','IX'],vaq)
            qiskit_circ.rx(2*error_rates_dict['Gcnot','XI'],aq+1)
            qiskit_circ.rxx(2*error_rates_dict['Gcnot','XX'],vaq,aq+1)            
    return qiskit_circ

def build_noisy_surface_code_qiskit(pcircuit, used_qubit_labels, qubit_dictionary,error_dict, incl_init=True):
    qiskit_circ=QuantumCircuit(len(used_qubit_labels))
    if incl_init:
        qiskit_circ.h(2)
        qiskit_circ.h(4)
        qiskit_circ.h(8)
        qiskit_circ.h(13)
        qiskit_circ.cx(8,3)
        qiskit_circ.cx(4,6)
        qiskit_circ.cx(13,15)
        qiskit_circ.cx(2,0)
        qiskit_circ.cx(8,6)
        qiskit_circ.cx(4,11)
        qiskit_circ.cx(3,2)
        qiskit_circ.cx(11,13)

    #syndrome extraction
    #use the pygsti circuit to build an expanded circuit w/ the correct number of rounds
    for idx,layer in enumerate(pcircuit):
        for lbl in layer:
            #print(lbl)
            if lbl[0] == 'Gh':
                # Factors of 2 b/c of difference between error gen rates and rotation angles
                qiskit_circ.h(qubit_dictionary[lbl[1]])
                qiskit_circ.unitary(error_dict['Gh'],[qubit_dictionary[lbl[1]]])
                ###Error rates used to be multiplied by 2. I'm editing it to see if it fixes discrepancies, but likely a bug elsewhere
                #qiskit_circ.rz(2*error_dict['Gh','Z'],qubit_dictionary[lbl[1]])
            elif lbl[0] == 'Gcnot':
                qiskit_circ.cx(qubit_dictionary[lbl[1]],qubit_dictionary[lbl[2]])
                qiskit_circ.unitary(error_dict['Gcnot'],[qubit_dictionary[lbl[2]], qubit_dictionary[lbl[1]]])
                #qiskit_circ.rx(2*error_dict['Gcnot','IX'],qubit_dictionary[lbl[2]])
                #qiskit_circ.rz(2*error_dict['Gcnot','ZI'],qubit_dictionary[lbl[1]])

    return qiskit_circ
    


def build_noisy_surface_code_qiskit_old(pcircuit, used_qubit_labels, qubit_dictionary,error_dict):
    #NOTE: only three kinds of error are used in these circuits right now

    # The circuit from https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.5.043137
    # initialization of logical 0 state
    qiskit_circ=QuantumCircuit(len(used_qubit_labels))
    qiskit_circ.h(2)
    qiskit_circ.h(4)
    qiskit_circ.h(8)
    qiskit_circ.h(13)
    qiskit_circ.cx(8,3)
    qiskit_circ.cx(4,6)
    qiskit_circ.cx(13,15)
    qiskit_circ.cx(2,0)
    qiskit_circ.cx(8,6)
    qiskit_circ.cx(4,11)
    qiskit_circ.cx(3,2)
    qiskit_circ.cx(11,13)

    #syndrome extraction
    #use the pygsti circuit to build an expanded circuit w/ the correct number of rounds
    for idx,layer in enumerate(pcircuit):
        for lbl in layer:
            #print(lbl)
            if lbl[0] == 'Gh':
                # Factors of 2 b/c of difference between error gen rates and rotation angles
                qiskit_circ.h(qubit_dictionary[lbl[1]])
                ###Error rates used to be multiplied by 2. I'm editing it to see if it fixes discrepancies, but likely a bug elsewhere
                qiskit_circ.rz(2*error_dict['Gh','Z'],qubit_dictionary[lbl[1]])
                qiskit_circ.rx(2*error_dict['Gh','X'],qubit_dictionary[lbl[1]])
            if lbl[0] == 'Gi':
                # Factors of 2 b/c of difference between error gen rates and rotation angles
                ###Error rates used to be multiplied by 2. I'm editing it to see if it fixes discrepancies, but likely a bug elsewhere
                qiskit_circ.rz(2*error_dict['Gi','Z'],qubit_dictionary[lbl[1]])
                qiskit_circ.rx(2*error_dict['Gi','X'],qubit_dictionary[lbl[1]])
            elif lbl[0] == 'Gcnot':
                qiskit_circ.cx(qubit_dictionary[lbl[1]],qubit_dictionary[lbl[2]])
                qiskit_circ.rx(2*error_dict['Gcnot','IX'],qubit_dictionary[lbl[2]])
                qiskit_circ.rz(2*error_dict['Gcnot','ZI'],qubit_dictionary[lbl[1]])
                qiskit_circ.rxx(2*error_dict['Gcnot','XX'],qubit_dictionary[lbl[2]], qubit_dictionary[lbl[1]])
                qiskit_circ.rzz(2*error_dict['Gcnot','ZZ'],qubit_dictionary[lbl[2]], qubit_dictionary[lbl[1]])

    return qiskit_circ

def sample_rep_code_detector_stats(d, sampler_circ, rounds, n_shots=10000000):
    sampler = StatevectorSampler()
    job = sampler.run([sampler_circ], shots=n_shots)
     
    # Extract the result for the 0th pub (this example only has one pub).
    result = job.result()[0]
    counts = result.data['meas'].get_counts()
    
    detector_counts = {}
    for bs, count in counts.items():
        #process the readout 
        readout = bs[::-1]
        ##THIS IS WRITTEN FOR THE REPETITION CODE##
        data_readout = np.array(int(readout[2*i]) for i in range(d))
        aux_readout = np.array([int(readout[2*i+1]) for i in range(d-1)]+[int(readout[i]) for i in range(2*d-1,len(readout))])
        #xor results of rounds together
        #first round is just as-is
        n_aux = d-1
        detectors = aux_readout[:n_aux]
        for j in range(1,rounds):
            detectors = np.append(detectors, np.bitwise_xor(aux_readout[n_aux*(j-1):n_aux*(j)], aux_readout[n_aux*(j):n_aux*(j+1)]))
        detector_string = ''.join([str(i) for i in detectors])
        if detector_string in detector_counts.keys():
            detector_counts[detector_string] += count
        else:
            detector_counts[detector_string] = count
    return detector_counts

def sample_surface_code_detector_stats(d, sampler_circ, aux_qs, data_qs, qubit_dictionary, rounds, n_shots=10000000):
    sampler = StatevectorSampler()
    job = sampler.run([sampler_circ], shots=n_shots)
     
    # Extract the result for the 0th pub (this example only has one pub).
    result = job.result()[0]
    counts = result.data['meas'].get_counts()
    
    detector_counts = {}
    for bs, count in counts.items():
        #process the readout 
        readout = bs[::-1]
        ##THIS IS WRITTEN FOR THE REPETITION CODE##
        #need to use the labels of the qubits!
        data_readout = np.array([int(readout[qubit_dictionary[j]]) for j in data_qs])
        aux_readout = np.array([int(readout[qubit_dictionary[j]]) for j in aux_qs])
        n_aux = d**2-1
        
        ##############
        #xor results of rounds together
        #first round is just as-is
        
        detectors = aux_readout[:n_aux]
        for j in range(rounds-1):
            detectors = np.append(detectors, np.bitwise_xor(aux_readout[n_aux*j:n_aux*(j+1)], aux_readout[n_aux*(j+1):n_aux*(j+2)]))
        detector_string = ''.join([str(i) for i in detectors])
        if detector_string in detector_counts.keys():
            detector_counts[detector_string] += count
        else:
            detector_counts[detector_string] = count
    return detector_counts

def compute_logical_z(state, mr_ancilla_qubits, data_qubits, qubit_dictionary):
    #get probability distribution for the ancilla qubits and the data qubit readout
    anc_qs = [qubit_dictionary[j] for j in mr_ancilla_qubits]
    data_qs = [qubit_dictionary[j] for j in data_qubits][:3]
    #print(mr_ancilla_qubits)
    # qargs = data_qs+anc_qs
    # print(qargs)
    # qargs.reverse()
    # print(qargs) 
    true_probs_readout = state.probabilities_dict()
    n_data = len(data_qubits)
    prob_dict_with_logical = Counter({})
    for readout, prob in true_probs_readout.items():
        logical_bits = [int(readout[-1*(d+1)]) for d in data_qs] ##hack
        logical = sum(logical_bits)%2 #TODO
        #print(readout, prob, logical_bits, logical)
        aux_readout = ''.join([str(readout[-1*(d+1)]) for d in anc_qs]) #Temporary hack that should work for d=1
        readout_w_logical = ''.join(aux_readout)+str(logical)
        prob_dict_with_logical.update({readout_w_logical: prob})
        #TODO sum everything with the same logical measurement result and same ancilla readout
        #add to new dictionary logical bit can be the last bit in the string
    return prob_dict_with_logical
    