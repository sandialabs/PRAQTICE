import numpy as np
import stim
import pygsti

import pygsti.tools.errgenproptools as eprop
from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator
from pygsti.baseobjs import Label, QubitSpace
from pygsti.models import LocalNoiseModel
from pygsti.modelmembers.operations import ComposedOp, LindbladErrorgen, ExpErrorgenOp, StaticCliffordOp

from pygsti.extras.dem_construction.pygsti_object_builders import *
from pygsti.extras.dem_construction.helper_functions import *

def setup_surface_code(d, rounds):
    
    stimcirc = stim.Circuit.generated("surface_code:rotated_memory_z", rounds=3, distance=d) #this is a dummy circuit
    num_qubits = stimcirc.num_qubits
    qubit_labels = range(num_qubits)
    
    pyg_c_temp = create_syndrome_extraction_circuit(stimcirc, qubit_labels, allowed_gates=['H', 'CX']) #looks like we don't know how to handle missing qubits
    
    x_ancilla_qubits = []
    coords = {}
    for line in stimcirc:
        if line.name == 'QUBIT_COORDS':
            coords[line.targets_copy()[0].value] = line.gate_args_copy()
        if line.name == 'H' and len(x_ancilla_qubits) == 0:
            x_ancilla_qubits += [target.value for target in line.targets_copy()]
        if line.name == 'MR':
            ancilla_qubits = [target.value for target in line.targets_copy()]
        if line.name == 'REPEAT':
            break
    used_qubit_labels = [int(q.split(') ')[1]) for q in str(stimcirc).split('R ')[0].split('\n')[:-1]]
    data_qubits = list(set(used_qubit_labels) - set(ancilla_qubits))
    z_ancilla_qubits = list(set(ancilla_qubits) - set(x_ancilla_qubits))
    extra_qs = list(set(qubit_labels)-set(used_qubit_labels))
    locations = [coords[bit] for bit in data_qubits]
    
    num_used_qubits = len(used_qubit_labels)

    # Builds all the pyGSTi objects we'll need to do the simulations.
    n_qubits = len(qubit_labels)
    n_ancilla = len(ancilla_qubits)
    mr_qubit_labels = [i for i in range(n_qubits+n_ancilla*(rounds-1))]
    mr_ancilla_qubits = ancilla_qubits+[i for i in range(n_qubits, n_qubits+n_ancilla*(rounds-1) )]

    mr_pcircuit = create_multiround_syndrome_extraction_circuit(stimcirc, qubit_labels, ancilla_qubits, rounds=rounds)

    return mr_pcircuit, mr_qubit_labels, mr_ancilla_qubits, data_qubits, ancilla_qubits, x_ancilla_qubits, z_ancilla_qubits, qubit_labels

def compute_tableau_and_inverse_tableau(mr_pcircuit, mr_qubit_labels, mr_pspec, data_qubits, z_ancilla_qubits, x_ancilla_qubits, qubit_labels, ancilla_qubits):
    z_stabilizer_qubits = []
    z_stabilizers = []
    num_expanded_qs = len(mr_qubit_labels)
    for q in z_ancilla_qubits:
        z_stabilizer_qubits.append([])
        for qs in mr_pspec.availability['Gcnot']:
            if q in qs:
                if qs[0] != q:
                    dq = qs[0]
                else:
                    dq = qs[1]
                z_stabilizer_qubits[-1].append(dq)
                
        ps = ['_']*num_expanded_qs
        for q in z_stabilizer_qubits[-1]:
            qi = mr_qubit_labels.index(q)
            ps[qi] = 'Z'
        ps = ''.join(ps)
        z_stabilizers.append(stim.PauliString(ps))
            
    x_stabilizer_qubits = []
    x_stabilizers = []
    for q in x_ancilla_qubits:
        x_stabilizer_qubits.append([])
        for qs in mr_pspec.availability['Gcnot']:
            if q in qs:
                if qs[0] != q:
                    dq = qs[0]
                else:
                    dq = qs[1]
                x_stabilizer_qubits[-1].append(dq)
        ps = ['_']*num_expanded_qs
        for q in x_stabilizer_qubits[-1]:
            qi = mr_qubit_labels.index(q)
            ps[qi] = 'X'
        ps = ''.join(ps)
        x_stabilizers.append(stim.PauliString(ps))
    
    logical_0_op = {}
    ps = ['_']*num_expanded_qs
    for q in data_qubits:
        qi = mr_qubit_labels.index(q)
        ps[qi] = 'Z'
    logical_0_op = stim.PauliString(ps)
    #logical_0_op = {}
    
    z_on_ancilla = []
    # Also includes Z on dummy qubits
    for q in mr_qubit_labels:
        if q not in data_qubits:
            z_on_ancilla.append(z_pauli_on_qubit(q, mr_qubit_labels))
            
    logical_0_stabilizers = z_stabilizers + x_stabilizers + [logical_0_op,] + z_on_ancilla

    output_state_stabilizers = []
    initial_state_tableau = []
    output_state_tableau = []
    output_state_inverse_tableau = []
    circuit_tableau = []
    
    sc = pygsti_c_to_stim(mr_pcircuit)
    circuit_tableau = sc.to_tableau(ignore_reset=True, ignore_measurement=True)
    initial_state_tableau = circuit_tableau.from_stabilizers(logical_0_stabilizers)
    for ps in logical_0_stabilizers:
        output_state_stabilizers.append(circuit_tableau(ps))
    output_state_tableau = circuit_tableau * initial_state_tableau
    output_state_inverse_tableau = output_state_tableau.inverse()
    
    #Check that Z is a stabilizer of single data qubits
    for q in data_qubits:
        assert(not is_stabilizer_or_antistabilizer_from_generators(z_pauli_on_qubit(q, mr_qubit_labels), output_state_stabilizers))  
    #Check that Z is a stabilizer of all the ancilla qubits ...
    for q in ancilla_qubits:
        assert(is_stabilizer_or_antistabilizer_from_generators(z_pauli_on_qubit(q, mr_qubit_labels), output_state_stabilizers))

    return output_state_tableau, output_state_inverse_tableau

def compute_stabilizers_data(mr_pcircuit, mr_qubit_labels, mr_pspec, data_qubits, z_ancilla_qubits, x_ancilla_qubits, qubit_labels, ancilla_qubits):
    z_stabilizer_qubits = []
    z_stabilizers = []
    num_expanded_qs = len(mr_qubit_labels)
    for q in z_ancilla_qubits:
        z_stabilizer_qubits.append([])
        for qs in mr_pspec.availability['Gcnot']:
            if q in qs:
                if qs[0] != q:
                    dq = qs[0]
                else:
                    dq = qs[1]
                z_stabilizer_qubits[-1].append(dq)
                
        ps = ['_']*num_expanded_qs
        for q in z_stabilizer_qubits[-1]:
            qi = mr_qubit_labels.index(q)
            ps[qi] = 'Z'
        ps = ''.join([p for i, p in enumerate(ps) if i in data_qubits])
        z_stabilizers.append(ps)
            
    x_stabilizer_qubits = []
    x_stabilizers = []
    for q in x_ancilla_qubits:
        x_stabilizer_qubits.append([])
        for qs in mr_pspec.availability['Gcnot']:
            if q in qs:
                if qs[0] != q:
                    dq = qs[0]
                else:
                    dq = qs[1]
                x_stabilizer_qubits[-1].append(dq)
        ps = ['_']*num_expanded_qs
        for q in x_stabilizer_qubits[-1]:
            qi = mr_qubit_labels.index(q)
            ps[qi] = 'X'
        ps = ''.join([p for i, p in enumerate(ps) if i in data_qubits])
        x_stabilizers.append(ps)

    return z_stabilizers, x_stabilizers
