'''
Utilities for model generation for 1Q-GST.

All functions are designed to be called from one rank.
'''

import inspect
import json
from os import stat
import numpy as np
import scipy as sp

import pygsti
from pygsti.modelmembers import operations as ops
from pygsti.modelmembers import states, povms


#Notes from discussions with KMR & SS on 2021-05-19, etc.
#1. strict over- or under-rotations (do nothing for idle, except perhaps a
#   z-rotation), selected from random.randn. May be maximally sensitive
#   to bare germ sets, see 1/L scaling for accuracy. Compare this to
#   depolarization operation noise.
#2. are there models that would violate the Heisenberg scaling, e.g., axis
#   misalignment error: pi/2*(cos(theta)x + sin(theta)y)?

#05-24
#3. separate parameters array/list into individual parameters for clarity.
#4. use argparse to input parameters into one_q_combined_gst.py.

#05-26
#5. for one-qubit systems, there are two possibe types of unitary errors: *axis*
#   of rotation (misalignment) and the *angle* of rotation
#6. consider the role of rotation angle and the overlap of rotation axes, along
#   with spam axes.

def create_noisy_model_from_lindblad_dict(target_model, errgen_dict, verbose=False):
    '''

    Parameters
    ----------
    

    Returns
    -------

    '''
    noisy_model = target_model.copy()
    for op_key, ideal_op in noisy_model.operations.items():
        if op_key == tuple() and '[]' in errgen_dict:
            # Hardcode empty tuple to "[]" since cannot JSON serialize tuples
            noise_dict = errgen_dict['[]']
        elif op_key in errgen_dict:
            noise_dict = errgen_dict[op_key]
        elif op_key.name in errgen_dict:
            noise_dict = errgen_dict[op_key.name]
        else:
            continue

        assert noise_dict is not None, "Assume have valid lindblad_term_dict by now"

        # TODO: We will have to be careful with this moving to 2-qubit,
        # but in theory this will let us do crosstalk noise as well
        errgen = ops.LindbladErrorgen.from_elementary_errorgens(noise_dict, state_space=noisy_model.state_space)
        experrop = ops.ExpErrorgenOp(errgen)
        noisy_op = ops.ComposedOp([ideal_op, experrop])

        noisy_model.operations[op_key] = noisy_op
    
    # SPAM error
    if "prep" in errgen_dict:
        for prep_key, ideal_prep in noisy_model.preps.items():
            noise_dict = errgen_dict["prep"]
            errgen = ops.LindbladErrorgen.from_elementary_errorgens(noise_dict, state_space=noisy_model.state_space)
            experrop = ops.ExpErrorgenOp(errgen)
            noisy_prep = states.ComposedState(ideal_prep.copy(), experrop)

            noisy_model.preps[prep_key] = noisy_prep
    
    if "povm" in errgen_dict:
        for povm_key, ideal_povm in noisy_model.povms.items():
            noise_dict = errgen_dict["povm"]
            errgen = ops.LindbladErrorgen.from_elementary_errorgens(noise_dict, state_space=noisy_model.state_space)
            experrop = ops.ExpErrorgenOp(errgen)
            noisy_povm = povms.ComposedPOVM(experrop, povm=ideal_povm.copy())

            noisy_model.povms[povm_key] = noisy_povm
    
    if verbose:
        print(noisy_model)

    return noisy_model

def combine_dicts(*args):
    '''

    Parameters
    ----------
    

    Returns
    -------

    '''
    combined = {}
    for d in args:
        for key, subdict in d.items():
            combined_subdict = combined.get(key, {})
            for subkey, subvalue in subdict.items():
                value = combined_subdict.get(subkey, 0)
                combined_subdict[subkey] = subvalue + value
            combined[key] = combined_subdict
    return combined

def sample_errgen_json(jsonfile, rand_state=None, verbose=False, num_qubits=1):
    noise_dict = {}
    linked_params = {}

    if rand_state is None or not isinstance(rand_state, np.random.RandomState):
        rand_state = np.random.RandomState(rand_state)

    with open(jsonfile, 'r') as f:
        errgen_spec = json.load(f)
        
    print(errgen_spec)
    
    for gatekey, spec_dict in errgen_spec.items():
        if '_' in gatekey:
            entries = gatekey.split('_')
            gatekey = (entries[0], int(entries[1]))
        noise_dict[gatekey] = {}

        # Expand "all H/S" macros
        basis = pygsti.baseobjs.Basis.cast('pp', 4**num_qubits)
        labels = basis.labels[1:] # Skip I
        expanded_spec_dict = {}
        for term, term_dict in spec_dict.items():
            if term == 'all H':
                if verbose: print('Expanding all H')
                for term_label in labels:
                    expanded_spec_dict[f'H{term_label}'] = term_dict
            elif term == 'all S':
                if verbose: print('Expanding all S')
                for term_label in labels:
                    expanded_spec_dict[f'S{term_label}'] = term_dict
            else:
                expanded_spec_dict[term] = term_dict

        for term, term_dict in expanded_spec_dict.items():
            if verbose:
                print('Sampling gate {}, term {}'.format(gatekey, term))

            try:
                sample_type = term_dict['type']
            except KeyError:
                if verbose:
                    print('  WARNING: Need to specify a "type" in term_dict, skipping')
                continue
            
            if sample_type == "normal":
                # Normal sampling with error scale
                noise_dict[gatekey][term] = term_dict["scale"]*rand_state.randn()
                if verbose:
                    print('  Sampled normal with error scale {} and got {}'.format(term_dict["scale"], noise_dict[gatekey][term]))
            elif sample_type == "uniform":
                interval = term_dict["interval"]
                # Sanity checks, only uniform bounds should make it here
                assert len(interval) == 2, "Expected a [low,high] bound for uniform sampling, got {}".format(interval)
                
                noise_dict[gatekey][term] = rand_state.uniform(interval[0], interval[1])
                if verbose:
                    print('  Sampled uniform with interval {} and got {}'.format(interval, noise_dict[gatekey][term]))
            elif sample_type == "static":
                noise_dict[gatekey][term] = term_dict["value"]
                if verbose:
                    print('  Set static value {}'.format(noise_dict[gatekey][term]))
            elif sample_type == "link":
                term_to_link = (term_dict["gate"], term_dict["term"])
                linked_params[gatekey, term] = term_to_link
                if verbose:
                    print('  Linking parameter for term {} with (gate, term) {}'.format(term, term_to_link))
            else:
                if verbose:
                    print('  WARNING: Unknown "type": {}, skipping'.format(sample_type))
    
    # Handle linked params once all non-linked params are sampled
    for gt1, gt2 in linked_params.items():
        g1, t1 = gt1
        g2, t2 = gt2
        assert t2 in noise_dict[g2], "{} not available for linking for {}".format(gt2, gt1)
        noise_dict[g1][t1] = noise_dict[g2][t2]
        if verbose:
            print('  Gate term {} linked with gate term {} with value {}'.format(gt1, gt2, noise_dict[g1][t1]))
    
    return noise_dict
