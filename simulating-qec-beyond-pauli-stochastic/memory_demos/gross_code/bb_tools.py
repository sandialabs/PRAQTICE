import stim
import numpy as np

from stimbposd import BPOSD, SinterDecoder_BPOSD, sinter_decoders

import time

from pygsti.extras.dem_construction import pygsti_object_builders as obj
from pygsti.extras.dem_construction import dem_tools as dems

import pygsti.tools.errgenproptools as eprop
from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator


l, m = 12, 6
N = l * m  # 72

def perm_matrix(dx, dy):
    P = np.zeros((N, N), dtype=np.uint8)
    for b in range(m):
        for a in range(l):
            i = gross_index(a, b)
            j = gross_index(a + dx, b + dy)
            P[j, i] = 1
    return P



def hstack_mod2(*mats): return np.concatenate(mats, axis=1) % 2
def vstack_mod2(*mats): return np.concatenate(mats, axis=0) % 2


# ---------- GF(2) math stuff ----------
def rref_mod2(m):
    m = m.copy() % 2
    n_rows, n_cols = m.shape
    pivots = []
    r = 0
    for c in range(n_cols):
        pivot = None
        for rr in range(r, n_rows):
            if m[rr, c]:
                pivot = rr; break
        if pivot is None: continue
        if pivot != r:
            m[[r, pivot]] = m[[pivot, r]]
        pivots.append(c)
        for rr in range(n_rows):
            if rr != r and m[rr, c]:
                m[rr, :] ^= m[r, :]
        r += 1
        if r == n_rows: break
    return m, pivots

def rank_mod2(m): return len(rref_mod2(m)[1])

def nullspace_basis_rows(m):
    m = m.copy() % 2
    R, piv = rref_mod2(m)
    n_rows, n_cols = R.shape
    piv_set = set(piv)
    free = [c for c in range(n_cols) if c not in piv_set]
    piv_to_row = {piv[i]: i for i in range(len(piv))}
    basis = []
    for fc in free:
        x = np.zeros(n_cols, dtype=np.uint8)
        x[fc] = 1
        for p in reversed(piv):
            row = piv_to_row[p]
            s = int(np.bitwise_and(R[row], x).sum() % 2)
            s_no_p = s ^ (R[row, p] & x[p])
            x[p] = s_no_p
        basis.append(x)
    return np.array(basis, dtype=np.uint8)

def rowspace_basis(m):
    B = np.zeros((0, m.shape[1]), dtype=np.uint8)
    for row in m:
        test = np.vstack([B, row[None, :]]) % 2
        if rank_mod2(test) > rank_mod2(B):
            B = test
    return B

def independent_mod_span(basis_rows, cand_row):
    if basis_rows.shape[0] == 0: return True
    old = rank_mod2(basis_rows)
    new = rank_mod2(np.vstack([basis_rows, cand_row[None, :]]))
    return new > old

# ---------- compute a logical-Z basis ----------
def compute_gross_logical_Z_basis():
    # check the ranks of the check operators
    # make sure there are 12 logicals
    rHX, rHZ = rank_mod2(HX), rank_mod2(HZ)
    assert rHX == 66 and rHZ == 66, (rHX, rHZ)
    k = 144 - rHX - rHZ
    assert k == 12

    # nullspace of HX
    K = nullspace_basis_rows(HX)           # 78 x 144
    S = rowspace_basis(HZ)                 # 66 x 144

    # pick L-only vectors first, then fill with mixed ones
    L_only = [v for v in K if v[N:].sum() == 0]
    mixed  = [v for v in K if v[N:].sum() != 0]

    basis = S.copy()
    picked = []

    for v in L_only:
        if independent_mod_span(basis, v):
            picked.append(v); basis = np.vstack([basis, v[None, :]])
        if len(picked) == 6: break

    for v in mixed:
        if independent_mod_span(basis, v):
            picked.append(v); basis = np.vstack([basis, v[None, :]])
        if len(picked) == 12: break

    assert len(picked) == 12 # Please I hope there are actually 12!
    return picked

# Convert those 12 supports into data-qubit indices consistent with the stim code above (L -> [72..143], R -> [144..215])
def supports_as_stim_indices(vectors):
    out = []
    for v in vectors:
        inds = []
        # L data
        for i in range(N):
            if v[i]: inds.append(N + i)
        # R data
        for i in range(N):
            if v[N + i]: inds.append(2*N + i)
        out.append(sorted(inds))
    return out

'''
Mostly follows Bravyi et al., High-threshold and low-overhead fault-tolerant quantum memory, 2024
'''

import stim
import numpy as np

def gross_index(a, b, l=12, m=6):
    return (a % l) + l * (b % m)

# Permutations for A and B
def perm_A1(i, l=12, m=6): b, a = divmod(i, l); return gross_index(a+3, b, l, m)
def perm_A2(i, l=12, m=6): b, a = divmod(i, l); return gross_index(a, b+1, l, m)
def perm_A3(i, l=12, m=6): b, a = divmod(i, l); return gross_index(a, b+2, l, m)
def perm_B1(i, l=12, m=6): b, a = divmod(i, l); return gross_index(a, b+3, l, m)
def perm_B2(i, l=12, m=6): b, a = divmod(i, l); return gross_index(a+1, b, l, m)
def perm_B3(i, l=12, m=6): b, a = divmod(i, l); return gross_index(a+2, b, l, m)

def flatten(nested):
    return [item for sublist in nested for item in sublist]

def build_gross_bb_repeated(rounds=3, l=12, m=6, p_cnot=None, p_init=None, p_meas=None, include_x_detectors=False):
    """
    Build a stim.Circuit for repeated syndrome extraction for the [[144,12,12]] bivariate bicycle code.
    """
    N = l * m
    qX, qL, qR, qZ = range(0, N), range(N, 2*N), range(2*N, 3*N), range(3*N, 4*N)

    full_c = stim.Circuit()

    def layer_cnot(pairs,circ):
        qs = flatten(pairs)
        circ.append("CNOT", qs)
        if p_cnot is not None:
            circ.append("DEPOLARIZE2", qs, p_cnot)
        #if p_cnot:
        #    for ctl, tgt in pairs:
        #        c.append("DEPOLARIZE2", [ctl, tgt], p_cnot)
        circ.append("TICK")
        return circ

    def pairs_from_perm(perm_fn, src_reg, dst_reg):
        return [(src_reg[perm_fn(i, l, m)], dst_reg[i]) for i in range(N)]

    prev_z_meas = None
    prev_x_meas = None

    ########

        # Init ancillas for this round
    full_c = stim.Circuit()
    full_c.append("RZ", qZ)
    if p_init:
        full_c.append("X_ERROR", qZ, p_init)
        #for q in qZ:
        #    c.append("X_ERROR", [q], p_init)
    full_c.append("RZ", qX)
    full_c.append('H',qX)
    if p_init:
        full_c.append("Z_ERROR", qX, p_init)
        #for q in qX:
        #    c.append("Z_ERROR", [q], p_init)
    full_c.append("TICK")

    # 7-layer CNOT schedule
    full_c = layer_cnot(pairs_from_perm(perm_A1, qR, qZ),full_c)
    full_c = layer_cnot(pairs_from_perm(perm_A2, qX, qL) + pairs_from_perm(perm_A3, qR, qZ),full_c)
    full_c = layer_cnot(pairs_from_perm(perm_B2, qX, qR) + pairs_from_perm(perm_B1, qL, qZ),full_c)
    full_c = layer_cnot(pairs_from_perm(perm_B1, qX, qR) + pairs_from_perm(perm_B2, qL, qZ),full_c)
    full_c = layer_cnot(pairs_from_perm(perm_B3, qX, qR) + pairs_from_perm(perm_B3, qL, qZ),full_c)
    full_c = layer_cnot(pairs_from_perm(perm_A1, qX, qL) + pairs_from_perm(perm_A2, qR, qZ),full_c)
    full_c = layer_cnot(pairs_from_perm(perm_A3, qX, qL),full_c)

    # Measure Z-checks
    if p_meas is not None:
        full_c.append("X_ERROR", qZ, p_meas)
        #for q in qZ:
        #    c.append("X_ERROR", [q], p_meas)
            
    full_c.append("MZ", qZ)

    # Detectors for Z stabilizers (compare with previous round)
    if prev_z_meas is not None:
        for i in range(N):
            full_c.append("DETECTOR", [stim.target_rec(-(3*N - i)), stim.target_rec(-(N - i))])
    else:
        for i in range(N):
            full_c.append("DETECTOR", [stim.target_rec(-(N - i))])
    prev_z_meas = [stim.target_rec(-(N - i)) for i in range(N)]
    full_c.append("RZ", qZ)
    if p_init is not None:
        full_c.append("X_ERROR", qZ, p_meas)
        #for q in qZ:
        #    c.append("X_ERROR", [q], p_meas)

    full_c.append("TICK")
    

    # Measure X-checks
    if p_meas is not None:
        full_c.append("Z_ERROR", qX, p_meas)
        #for q in qX:
        #    c.append("Z_ERROR", [q], p_meas)
    full_c.append('H',qX)            
    full_c.append("MZ", qX)
    if prev_x_meas is not None:
        if include_x_detectors:
            for i in range(N):
                full_c.append("DETECTOR", [stim.target_rec(-(3*N - i)), stim.target_rec(-(N - i))])
                #c.append("DETECTOR", [prev_x_meas[i], stim.target_rec(-(N - i))])
    else:
        if include_x_detectors:
            for i in range(N):
                full_c.append("DETECTOR", [stim.target_rec(-(N - i))])
    prev_x_meas = [stim.target_rec(-(N - i)) for i in range(N)]
    full_c.append("RZ", qX)
    full_c.append('H',qX)
    if p_init is not None:
        full_c.append("Z_ERROR", qX, p_meas)
        #for q in qX:
        #    c.append("Z_ERROR", [q], p_meas)

    full_c.append("TICK")


    ########
        # Init ancillas for this round
    c = stim.Circuit()
    c.append("RZ", qZ)
    if p_init:
        c.append("X_ERROR", qZ, p_init)
        #for q in qZ:
        #    c.append("X_ERROR", [q], p_init)
    c.append("RZ", qX)
    c.append('H',qX)
    if p_init:
        c.append("Z_ERROR", qX, p_init)
        #for q in qX:
        #    c.append("Z_ERROR", [q], p_init)
    c.append("TICK")

    # 7-layer CNOT schedule
    c = layer_cnot(pairs_from_perm(perm_A1, qR, qZ),c)
    c = layer_cnot(pairs_from_perm(perm_A2, qX, qL) + pairs_from_perm(perm_A3, qR, qZ),c)
    c = layer_cnot(pairs_from_perm(perm_B2, qX, qR) + pairs_from_perm(perm_B1, qL, qZ),c)
    c = layer_cnot(pairs_from_perm(perm_B1, qX, qR) + pairs_from_perm(perm_B2, qL, qZ),c)
    c = layer_cnot(pairs_from_perm(perm_B3, qX, qR) + pairs_from_perm(perm_B3, qL, qZ),c)
    c = layer_cnot(pairs_from_perm(perm_A1, qX, qL) + pairs_from_perm(perm_A2, qR, qZ),c)
    c = layer_cnot(pairs_from_perm(perm_A3, qX, qL),c)

    # Measure Z-checks
    if p_meas is not None:
        c.append("X_ERROR", qZ, p_meas)
        #for q in qZ:
        #    c.append("X_ERROR", [q], p_meas)
            
    c.append("MZ", qZ)

    # Detectors for Z stabilizers (compare with previous round)
    if prev_z_meas is not None:
        for i in range(N):
            c.append("DETECTOR", [stim.target_rec(-(3*N - i)), stim.target_rec(-(N - i))])
    else:
        for i in range(N):
            c.append("DETECTOR", [stim.target_rec(-(N - i))])
    prev_z_meas = [stim.target_rec(-(N - i)) for i in range(N)]
    c.append("RZ", qZ)
    if p_init is not None:
        c.append("X_ERROR", qZ, p_meas)
        #for q in qZ:
        #    c.append("X_ERROR", [q], p_meas)

    c.append("TICK")

    # Measure X-checks
    if p_meas is not None:
        c.append("Z_ERROR", qX, p_meas)
        #for q in qX:
        #    c.append("Z_ERROR", [q], p_meas)
    c.append("H", qX)        
    c.append("MZ", qX)
    if prev_x_meas is not None:
        if include_x_detectors:
            for i in range(N):
                c.append("DETECTOR", [stim.target_rec(-(3*N - i)), stim.target_rec(-(N - i))])
                #c.append("DETECTOR", [prev_x_meas[i], stim.target_rec(-(N - i))])
    else:
        if include_x_detectors:
            for i in range(N):
                c.append("DETECTOR", [stim.target_rec(-(N - i))])
    prev_x_meas = [stim.target_rec(-(N - i)) for i in range(N)]
    c.append("RZ", qX)
    c.append("H", qX)
    if p_init is not None:
        c.append("Z_ERROR", qX, p_meas)
        #for q in qX:
        #    c.append("Z_ERROR", [q], p_meas)

    c.append("TICK")

    ####
    for j in range(rounds-1):
        full_c.append(c)

    #add final measurements to the circuit
    basis = compute_gross_logical_Z_basis()
    stim_supports = supports_as_stim_indices(basis) #these are the qubit indices for the logical ops
    
    if p_meas is not None:
        for q in qZ:
            full_c.append("X_ERROR", [q], p_meas)
    full_c.append('MZ', [i for i in range(N,3*N)]) #measure the data qubits

    pA3 = pairs_from_perm(perm_A3, qR, qZ)
    pA1 = pairs_from_perm(perm_A1, qR, qZ)
    pB1 = pairs_from_perm(perm_B1, qL, qZ)
    pB2 = pairs_from_perm(perm_B2, qL, qZ)
    pB3 = pairs_from_perm(perm_B3, qL, qZ)
    pA2 = pairs_from_perm(perm_A2, qR, qZ)
    stabilizers = list([t[0] for t in k] for k in list(zip(pA3, pA1, pB1, pB2, pB3, pA2)))

    for j, indices in enumerate(stabilizers):
        targs = [stim.target_rec(-(3*N-i)) for i in indices]+[stim.target_rec(-(4*N-j))]
        full_c.append('DETECTOR', targs)

    for j, indices in enumerate(stim_supports):
        full_c.append('OBSERVABLE_INCLUDE', arg=[j], targets=[stim.target_rec(-(3*N-i)) for i in indices]) #subtract N because we're only measuring the data qubits
    
    return full_c

# A = x^3 + y + y^2 , B = y^3 + x + x^2
A = (perm_matrix(3, 0) ^ perm_matrix(0, 1) ^ perm_matrix(0, 2))  # 72 x 72
B = (perm_matrix(0, 3) ^ perm_matrix(1, 0) ^ perm_matrix(2, 0))  # 72 x 72

HX = hstack_mod2(A, B)            # 72 x 144
HZ = hstack_mod2(B.T, A.T)        # 72 x 144

basis = compute_gross_logical_Z_basis()
stim_supports = supports_as_stim_indices(basis)

