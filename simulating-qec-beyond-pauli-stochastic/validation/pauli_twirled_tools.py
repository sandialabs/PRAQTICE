import pygsti
import numpy as np
import stim 
from scipy.linalg import hadamard

def add_error_to_stim_circuit(circuit, stim_edict, measurement_names, reset_names, include_idles=True, include_meas_idles=False, init_circuit=None):
    '''
    adds error to a stim circuit
    stim_edict is a dictionary whose keys are stim gate names and whose values are tuples of the form (stim_error_name, rate).
    Here, stim_error_name is the stim name for the error channel, e.g. 'X_ERROR'. 
    include_idles allows you to pad layers of a circuit with idles (which you can then add error to)
    '''
    if init_circuit is None:
        noisy_c = stim.Circuit()
    else:
        noisy_c = init_circuit.copy()
    ###If gate, measurement, or reset: add error according to dictionary###
    for line in circuit:
        if type(line)== stim.CircuitRepeatBlock:
            new_body = stim.Circuit()
            for l in line.body_copy():
                op = l.name
                targs = l.targets_copy()
                if op in stim_edict.keys():
                    error_layer = stim_edict[op]
                    if op not in measurement_names and op not in reset_names:
                        new_body.append(l)
                        new_body.append(error_layer[0],targs, error_layer[1])
                        if include_idles:
                            unused_qs = [k for k in range(circuit.num_qubits) if k not in [j.value for j in l.targets_copy()]]
                            error_layer = stim_edict['I']
                            new_body.append(error_layer[0], unused_qs, error_layer[1])
                            #print(line.targets_copy(), unused_qs)
                    else: #note that 'MR' is a measurement and a reset
                        if op in measurement_names:
                            new_body.append(error_layer[0],targs, error_layer[1])
                            new_body.append(l)
                        if op in reset_names:
                            new_body.append(l)
                            new_body.append(error_layer[0],targs, error_layer[1])
                        if include_meas_idles:
                            unused_qs = [k for k in range(circuit.num_qubits) if k not in [j.value for j in l.targets_copy()]]
                            error_layer = stim_edict['I']
                            new_body.append(error_layer[0], unused_qs, error_layer[1])
                            #print(line.targets_copy(), unused_qs)

                else:
                    new_body.append(l)
                        
    
            noisy_c.append(stim.CircuitRepeatBlock(repeat_count=line.repeat_count, body=new_body))
        
        else:
            op = line.name
            targs = line.targets_copy()
            if op in stim_edict.keys():
                error_layer = stim_edict[op]
                if op not in measurement_names and op not in reset_names:
                    noisy_c.append(line)
                    noisy_c.append(error_layer[0],targs, error_layer[1])
                    if include_idles:
                        unused_qs = [k for k in range(circuit.num_qubits) if k not in [j.value for j in line.targets_copy()]]
                        error_layer = stim_edict['I']
                        noisy_c.append(error_layer[0], unused_qs, error_layer[1])
                        #print(line.targets_copy(), unused_qs)
                else: #note that 'MR' is a measurement and a reset
                    if op in measurement_names:
                        noisy_c.append(error_layer[0],targs, error_layer[1])
                        noisy_c.append(line)
                    if op in reset_names:
                        noisy_c.append(line)
                        noisy_c.append(error_layer[0],targs, error_layer[1])
                    if include_meas_idles:
                        unused_qs = [k for k in range(circuit.num_qubits) if k not in [j.value for j in line.targets_copy()]]
                        error_layer = stim_edict['I']
                        noisy_c.append(error_layer[0], unused_qs, error_layer[1])
                    

            else:
                noisy_c.append(line)
                    
    return noisy_c

def build_pauli_twirled_model(model, circuit, spam_error=0, include_idles=False, include_meas_idles=False, init_circuit=None, return_c=False):
    '''
    model: pyGSTi model for the noisy gates
    circuit: Stim circuit with no error, to add the error to. Should have the desired detector structure for the DEM.
    returns: stim DEM for noisy circuit
    Note to self: Have to be careful if we want to simulate a circuit that has ideal state prep followed by noisy SE rounds
    '''
    #CNOT#
    twoQ_rates_order  = ['II', 'IY', 'IX', 'IZ', 'YI', 'YY', 'YX', 'YZ', 'XI', 'XY', 'XX', 'XZ', 'ZI', 'ZY', 'ZX', 'ZZ']
    pauli_eigs_cnot = np.diag(model['Gcnot']) #model.operation_blks[pygsti.baseobjs.label.Label('Gcnot',(0,1))].factorops[1].to_dense())
    rates = hadamard(16) @ pauli_eigs_cnot /16 #II, IY, IX, IZ, YI , YY YX, YZ, XI, XY, XX, XZ, ZI, ZY, ZX, ZZ   format: P_Q1 P_Q2
    stim_order = ['IX', 'IY', 'IZ', 'XI', 'XX', 'XY', 'XZ', 'YI', 'YX', 'YY', 'YZ', 'ZI', 'ZX', 'ZY', 'ZZ']
    stim_rates_cnot = [rates[twoQ_rates_order.index(p)] if rates[twoQ_rates_order.index(p)]>1e-12 else 0 for p in stim_order] #IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ
    #H gate#
    pauli_eigs_h = np.diag(model['Gh']) #model.operation_blks[pygsti.baseobjs.label.Label('Gh',(0,))].factorops[1].to_dense())
    pauli_rates_h = hadamard(4) @ pauli_eigs_h /4 #I, Y, X, Z
    stim_rates_h = [pauli_rates_h[2], pauli_rates_h[1], pauli_rates_h[3]]
    stim_rates_h = [r if r>1e-12 else 0 for r in stim_rates_h]

    #I gate#
    pauli_eigs_i = np.diag(model['Gi']) #model.operation_blks[pygsti.baseobjs.label.Label('Gi',(0,))].factorops[1].to_dense())
    pauli_rates_i = hadamard(4) @ pauli_eigs_i /4 #I, Y, X, Z
    stim_rates_i = [pauli_rates_i[2], pauli_rates_i[1], pauli_rates_i[3]]
    stim_rates_i = [r if r>1e-12 else 0 for r in stim_rates_i]

    #stim error dictionary
    stim_edict = {}
    #the listed numbers are Pauli error probabilities for a post-gate error model. This is just an example. 
    stim_edict['CX']= ('PAULI_CHANNEL_2', stim_rates_cnot)
    stim_edict['H']= ('PAULI_CHANNEL_1', stim_rates_h)
    stim_edict['I']= ('PAULI_CHANNEL_1', stim_rates_i) #this seems low enough thtat it isn't worth scaling down

    stim_edict['M']= ('X_ERROR',spam_error)
    stim_edict['R']= ('X_ERROR',0)

    #stim circuit#
    stimc = add_error_to_stim_circuit(circuit, stim_edict, ['M',  'MR'], ['R'] ,include_idles=include_idles, include_meas_idles=include_meas_idles, init_circuit=init_circuit)
    #stim DEM#
    stim_dem = stimc.detector_error_model(approximate_disjoint_errors=True)
    if return_c:
        return stim_dem, stimc
    else:
        return stim_dem

def unroll_repeat_blocks(circuit, add_2q_idles=False, add_measurement_idles=False):
    new_circuit = stim.Circuit()
    
    for operation in circuit:
        if operation.name == "REPEAT":
            n_repeats = operation.repeat_count  # Assuming the first arg is the repeat count
            for _ in range(n_repeats):
                # Append the operations that follow the repeat block
                for op in operation.body_copy():
                    if op.name == "REPEAT":
                        break  # Stop if we hit another repeat block
                    new_circuit.append_operation(op.name, op.targets_copy(), op.gate_args_copy())
                    if (op.name=='MR' or op.name=='R') and add_measurement_idles:
                        #add idles on unmeasured qubits
                        unmeasured_qubits = set([q for q in range(circuit.num_qubits)])-set(t.qubit_value for t in op.targets_copy())
                        #print(op, unmeasured_qubits)
                        new_circuit.append_operation('I', unmeasured_qubits, op.gate_args_copy())
                    elif (op.name=='CX' or op.name=='CZ') and add_2q_idles:
                        unused_qubits = set([q for q in range(circuit.num_qubits)])-set(t.qubit_value for t in op.targets_copy())
                        #print(op, unused_qubits)
                        new_circuit.append_operation('I', unused_qubits, op.gate_args_copy())
        else:
            # Append operations that are not part of a repeat block
            new_circuit.append_operation(operation.name, operation.targets_copy())
            if (operation.name=='MR' or operation.name=='R') and add_measurement_idles:
                unmeasured_qubits = set([q for q in range(circuit.num_qubits)])-set(t.qubit_value for t in operation.targets_copy())
                #print(operation, unmeasured_qubits)
                new_circuit.append_operation('I', unmeasured_qubits, operation.gate_args_copy())
            elif (operation.name=='CX' or operation.name=='CZ') and add_2q_idles:
                unused_qubits = set([q for q in range(circuit.num_qubits)])-set(t.qubit_value for t in operation.targets_copy())
                #print(operation, unused_qubits)
                new_circuit.append_operation('I', unused_qubits, operation.gate_args_copy())
    
    return new_circuit

