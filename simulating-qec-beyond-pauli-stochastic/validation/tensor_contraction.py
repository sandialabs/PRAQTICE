import stim
import quimb as qu
import quimb.tensor as qtn
from itertools import product
import io
import numpy as np

def xor_tensor(num_vars, target_bit):
    # Construct a tensor that represents the xor operator
    shape = (2,) * num_vars
    data = []
    for bits in product([0, 1], repeat=num_vars):
        val = sum(bits) % 2
        data.append(1.0 if val == target_bit else 0.0) 
    return np.array(data).reshape(shape)

def stim_dem_to_tensor_network(dem: stim.DetectorErrorModel, target_bits):
    assert len(target_bits) == dem.num_detectors, "Number of target bits must equal number of detectors"
    tn = qtn.TensorNetwork([])
    detector_to_errors = {}
    error_index_map = {}  
    

    # tensor for each DEM event
    for i, event in enumerate(dem):
        prob = event.args_copy()[0]
        targets = [t.val for t in event.targets_copy() if t.is_relative_detector_id]
        main_ind = f'e{i}_main'
        prob_tensor = qtn.Tensor(
            data=[1 - prob, prob],
            inds=(main_ind,),
            tags={f'p{i}'}
        )
        tn |= prob_tensor

        aux_inds = [f'e{i}_{j}' for j in range(len(targets))]
        error_index_map[i] = aux_inds

        all_inds = [main_ind] + aux_inds
        id_tensor = np.zeros((2,) * len(all_inds))
        for bits in product([0, 1], repeat=len(all_inds)):
            if all(b == bits[0] for b in bits):
                id_tensor[bits] = 1.0
        tn |= qtn.Tensor(data=id_tensor, inds=all_inds, tags={f'e{i}'})

        for idx, d in enumerate(targets):
            detector_to_errors.setdefault(d, []).append((i, aux_inds[idx]))

    # XOR constraint tensors for each detector
    for d, target in enumerate(target_bits):
        if d not in detector_to_errors:
            if target != 0:
                return 0.0  # impossible
            continue

        involved = detector_to_errors[d]
        inds = [ind for (_, ind) in involved]
        t = xor_tensor(len(inds), target)
        tn |= qtn.Tensor(data=t, inds=inds, tags={f'd{d}'})

    return tn


def compute_probability_from_dem(dem: stim.DetectorErrorModel, target_bits):
    tn = stim_dem_to_tensor_network(dem, target_bits)
    return tn, tn.contract(all, optimize='auto-hq')

def compute_all_probs(dem: stim.DetectorErrorModel):
    prob_dict = {}
    for bits in product([0, 1], repeat=dem.num_detectors):
        tn, prob = compute_probability_from_dem(dem, bits)
        prob_dict[''.join([str(j) for j in bits])] = prob
    return prob_dict