"""CUDA-Q implementation of QUOPS for the quops.ipynb tutorial.

Sections: constants; helper functions; circuit construction; randomized
compilation; MCFE caps; noise models; statistical inference; scores and rates;
scan scheduling; plotting; reproducible physical acquisition; logical
resource displays; inference calibration; command-line tools.

The notebook keeps its CUDA-Q demo definitions visible and explicitly imports
the supporting functions it uses. All regression tests live in tests.py.
Logical programs stay visible in the notebook; the helpers here display their
already-compiled resource counts without requiring Logical at module import.

QUOPS estimates mean process polarization over a random circuit ensemble.
Its mirror estimator relies on the reference-compiler assumptions discussed
in the paper. The implemented ratio of ensemble means is approximate;
sampling confidence bounds do not remove its model or normalization bias.
The tutorial uses the paper's gated score tests and Hochberg capability-region
tests with normal/bootstrap uncertainty. Legacy fixed-family archive analysis
retains its Bonferroni adjustment and conservative degenerate-data fallback.

Run a physical scan with ``python src.py --help``; run the known-channel
inference calibration with ``python src.py calibrate --help``.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
import platform
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterable, Sequence

import cudaq
import numpy as np
import pandas as pd


# ===========================================================================
# Constants and gate conventions
# ===========================================================================

RX, RY, RZ, X, Y, Z, CX, U3 = range(8)

ROTATION_GATES = (RX, RY, RZ)
TWO_QUBIT_GATES = (CX,)

GATE_NAMES = {RX: "rx", RY: "ry", RZ: "rz", X: "x", Y: "y", Z: "z", CX: "cx", U3: "u3"}


# Pauli indices used throughout: 0=I, 1=Z, 2=X, 3=Y, so that the low bit is the
# Z component and the high bit is the X component.
PAULI_I, PAULI_Z, PAULI_X, PAULI_Y = range(4)
PAULI_GATE_ID = (None, Z, X, Y)

_PAULI_MATRICES = (np.eye(2, dtype=complex), np.diag([1.0, -1.0]).astype(complex), np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex), np.array([[0.0, -1j], [1j, 0.0]], dtype=complex))


# ===========================================================================
# Helper functions: scalar validation and single-qubit algebra
# ===========================================================================

def _integer_value(value, name: str) -> int:
    """Coerce an integer-valued scalar without silently truncating it."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    if isinstance(value, (int, np.integer)):
        return int(value)
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{name} must be an integer")
    return int(numeric)

def _positive_integer(name: str, value) -> int:
    value = _integer_value(value, name)
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value

def _nonnegative_integer(name: str, value) -> int:
    value = _integer_value(value, name)
    if value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return value

def _probability(name: str, value: float) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be numeric, not boolean")
    value = float(value)
    if not math.isfinite(value) or not 0.0 < value < 1.0:
        raise ValueError(f"{name} must be finite and between 0 and 1")
    return value

def _noise_probability(name: str, value: float) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be numeric, not boolean")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite probability in [0, 1]") from exc
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be a finite probability in [0, 1]")
    return value


def pauli_matrix(index: int) -> np.ndarray:
    """Return the 2x2 matrix for a Pauli index (0=I, 1=Z, 2=X, 3=Y)."""
    return _PAULI_MATRICES[int(index)]


def conjugate_pauli_through_cnot(control: int, target: int) -> tuple[int, int]:
    """Return ``CX (P_c ⊗ P_t) CX`` as a pair of Pauli indices.

    Under conjugation by CNOT: ``X_c -> X_c X_t``, ``Z_t -> Z_c Z_t``, while
    ``Z_c`` and ``X_t`` are unchanged.  Because CNOT is Clifford, a Pauli frame
    stays a Pauli frame; this exact closure is what lets randomized compilation
    move compensating twirls through CNOTs without changing the ideal unitary.
    """
    z_c, x_c = control & 1, control >> 1
    z_t, x_t = target & 1, target >> 1
    out_control = x_c << 1 | (z_c ^ z_t)
    out_target = (x_c ^ x_t) << 1 | z_t
    return out_control, out_target


def rz_matrix(angle: float) -> np.ndarray:
    return np.diag([np.exp(-0.5j * angle), np.exp(0.5j * angle)])


def ry_matrix(angle: float) -> np.ndarray:
    cosine, sine = math.cos(angle / 2.0), math.sin(angle / 2.0)
    return np.array([[cosine, -sine], [sine, cosine]], dtype=complex)


def rx_matrix(angle: float) -> np.ndarray:
    cosine, sine = math.cos(angle / 2.0), math.sin(angle / 2.0)
    return np.array([[cosine, -1j * sine], [-1j * sine, cosine]], dtype=complex)


def u3_matrix(theta: float, phi: float, lam: float) -> np.ndarray:
    """Return ``Rz(phi) Ry(theta) Rz(lam)``, CUDA-Q's ``u3`` up to global phase."""
    return rz_matrix(phi) @ ry_matrix(theta) @ rz_matrix(lam)


def single_qubit_matrix(gate: int, angle: float, phi: float = 0.0, lam: float = 0.0):
    """Return the 2x2 matrix for any single-qubit gate id in this module."""
    if gate == RX:
        return rx_matrix(angle)
    if gate == RY:
        return ry_matrix(angle)
    if gate == RZ:
        return rz_matrix(angle)
    if gate == X:
        return pauli_matrix(PAULI_X)
    if gate == Y:
        return pauli_matrix(PAULI_Y)
    if gate == Z:
        return pauli_matrix(PAULI_Z)
    if gate == U3:
        return u3_matrix(angle, phi, lam)
    raise ValueError(f"not a single-qubit gate id: {gate}")


def u3_angles_from_matrix(unitary: np.ndarray) -> tuple[float, float, float]:
    """Decompose a 2x2 unitary into ``(theta, phi, lam)`` up to global phase.

    Inverse of :func:`u3_matrix`, so ``u3_matrix(*u3_angles_from_matrix(U))``
    equals ``U`` up to a global phase.
    """
    unitary = np.asarray(unitary, dtype=complex)
    a, b = unitary[0, 0], unitary[0, 1]
    c, d = unitary[1, 0], unitary[1, 1]

    # With theta in [0, pi] both cos(theta/2) and sin(theta/2) are non-negative,
    # so the entry magnitudes fix theta and the entry phases fix phi and lam.
    theta = 2.0 * math.atan2(abs(c), abs(a))

    # phi and lam are each recovered relative to a common entry, which keeps the
    # relative signs consistent.  Solving instead for phi + lam and phi - lam
    # separately loses a factor of -1 whenever either wraps past pi, because Rz
    # is 4 pi periodic while numpy.angle is 2 pi periodic.
    if abs(a) < 1e-12:
        # theta = pi: only phi - lam is determined, so pin lam to zero.
        return theta, float(np.angle(c) - np.angle(-b)), 0.0
    if abs(c) < 1e-12:
        # theta = 0: only phi + lam is determined, so pin lam to zero.
        return theta, float(np.angle(d) - np.angle(a)), 0.0

    phi = float(np.angle(c) - np.angle(a))
    lam = float(np.angle(-b) - np.angle(a))
    return theta, phi, lam


# ===========================================================================
# Circuit construction: CUDA-Q kernels and mirror composition
# ===========================================================================
#
# Implements the QUOPS circuit ensemble of Supplementary Section II B: a
# width-``w`` circuit consists of benchmark layers, each a pair containing a
# CNOT layer followed by an ``R_P(theta)`` layer.  The paper calls the total
# number of elementary gate layers ``d``, so it contains ``d / 2`` benchmark
# layers. All public ``depth`` arguments use the paper's even ``d``.
# Private *_from_layer_pairs helpers handle explicitly legacy archive units.
# The final benchmark layer is partially
# filled, acting on exactly ``zeta * w`` uniformly chosen qubits.
#
# Circuit size follows Supplementary Section I A: one per single-qubit gate and
# two per CNOT.  This is the architecture-independent logical size ``s`` of the
# sampled QUOPS circuit, not the number of native gates after routing or gate
# synthesis.  Compilation overhead is intentionally charged through the
# observed polarization and runtime of the implementation being benchmarked.
#
# Circuits are carried as flat gate arrays alongside the compiled CUDA-Q
# kernel.  CUDA-Q exposes no stable public API for recovering a gate list from
# a compiled kernel, and the QUOPS protocol needs to rewrite circuits
# (adjoints, randomized compilation, concatenation), so the arrays are the
# source of truth and the kernel is derived from them.
#
# Exact versus experimental objects:
#
# * The sampled arrays define the ideal circuit C exactly.
# * ``build_kernel`` realizes that ideal manifest in CUDA-Q.
# * A backend may compile C to a different native circuit g(C); QUOPS measures
#   the quality of that complete implementation, including compilation error.
# * The mirror circuits must remain physically present.  Optimizing C followed
#   by C^dag to the identity would erase the errors MCFE is meant to observe.


def build_kernel(gates, q0, q1, angles, width, phis=None, lams=None):
    """Build a CUDA-Q kernel from flat gate arrays and attach those arrays.

    ``U3(theta, phi, lam)`` uses ``angles`` for ``theta`` and the optional
    ``phis``/``lams`` arrays for the remaining Euler angles; every other gate
    ignores them.

    The arrays are copied before the CUDA-Q function is defined so that the
    values captured by the compiler stay identical to the arrays attached to
    the returned kernel even if a caller later mutates its own lists.

    ``atomic_quantum_region=True`` is scientifically important for MCFE, not
    merely an implementation detail.  M1 and M2 intentionally contain a
    forward circuit adjacent to an inverse.  Their invocation boundaries must remain
    intact so cross-component cancellation cannot erase benchmark errors.
    CUDA-Q can still optimize inside each component; atomic regions do not
    promise literal preservation of every internal gate.
    """
    # CUDA-Q permits runtime indexing into captured lists but requires tuple
    # indices to be compile-time constants, so private list copies give both
    # compiler compatibility and isolation from later caller mutation.
    width = _integer_value(width, "width")
    if width < 1:
        raise ValueError("width must be at least one")

    frozen_gates = [_integer_value(gate, "gate identifier") for gate in gates]
    frozen_q0 = [_integer_value(qubit, "q0 index") for qubit in q0]
    frozen_q1 = [_integer_value(qubit, "q1 index") for qubit in q1]
    frozen_angles = [float(angle) for angle in angles]
    n = len(frozen_gates)
    frozen_phis = [0.0] * n if phis is None else [float(value) for value in phis]
    frozen_lams = [0.0] * n if lams is None else [float(value) for value in lams]

    if not len(frozen_q0) == len(frozen_q1) == len(frozen_angles) == len(frozen_phis) == len(frozen_lams) == n:
        raise ValueError("gates, q0, q1, angles, phis, and lams must have the same length")

    unknown_gates = sorted(set(frozen_gates) - set(GATE_NAMES))
    if unknown_gates:
        raise ValueError(f"unknown gate identifiers: {unknown_gates}")
    if not all(math.isfinite(value) for value in (*frozen_angles, *frozen_phis, *frozen_lams)):
        raise ValueError("gate angles must be finite")
    for index, (gate, first, second) in enumerate(zip(frozen_gates, frozen_q0, frozen_q1)):
        if not 0 <= first < width:
            raise ValueError(f"q0[{index}]={first} is outside width {width}")
        if gate == CX:
            if not 0 <= second < width:
                raise ValueError(f"q1[{index}]={second} is outside width {width}")
            if first == second:
                raise ValueError(f"CX at index {index} must use distinct qubits")

    if n == 0:
        # CUDA-Q cannot infer the element type of an empty captured list even
        # when the loop trip count is zero, so represent the identity directly.
        @cudaq.kernel(atomic_quantum_region=True)
        def C(q: cudaq.qview):
            pass

    else:

        @cudaq.kernel(atomic_quantum_region=True)
        def C(q: cudaq.qview):
            for i in range(n):
                g = frozen_gates[i]
                if g == RX:
                    rx(frozen_angles[i], q[frozen_q0[i]])
                elif g == RY:
                    ry(frozen_angles[i], q[frozen_q0[i]])
                elif g == RZ:
                    rz(frozen_angles[i], q[frozen_q0[i]])
                elif g == X:
                    x(q[frozen_q0[i]])
                elif g == Y:
                    y(q[frozen_q0[i]])
                elif g == Z:
                    z(q[frozen_q0[i]])
                elif g == CX:
                    x.ctrl(q[frozen_q0[i]], q[frozen_q1[i]])
                elif g == U3:
                    u3(frozen_angles[i], frozen_phis[i], frozen_lams[i], q[frozen_q0[i]])


    C.gates = frozen_gates
    C.q0s = frozen_q0
    C.q1s = frozen_q1
    C.angles = frozen_angles
    C.phis = frozen_phis
    C.lams = frozen_lams
    C.width = width
    return C


def require_kernel_arrays(kernel):
    """Return ``(gates, q0s, q1s, angles, phis, lams, width)`` for a built kernel."""
    if not hasattr(kernel, "gates"):
        raise ValueError("expected a kernel built by build_kernel")
    return kernel.gates, kernel.q0s, kernel.q1s, kernel.angles, kernel.phis, kernel.lams, int(kernel.width)


def circuit_size(gates: Iterable[int]) -> int:
    """Return the logical QUOPS size: one per 1Q gate and two per CNOT.

    This weighting is the paper's architecture-neutral operation count.  It
    makes a fully populated benchmark layer have size approximately ``2 * w``:
    ``w`` from the rotation layer and ``2 * floor(w / 2)`` from its CNOTs.
    """
    gate_ids = [_integer_value(gate, "gate identifier") for gate in gates]
    unknown_gates = sorted(set(gate_ids) - set(GATE_NAMES))
    if unknown_gates:
        raise ValueError(f"unknown gate identifiers: {unknown_gates}")
    return sum(2 if gate in TWO_QUBIT_GATES else 1 for gate in gate_ids)


def adjoint_C(kernel):
    """Return the exact ideal inverse ``C^dag`` used by M1 and M2.

    Taking an adjoint reverses temporal order and inverts every gate.  CNOT and
    Pauli gates are self-inverse; axial rotations negate their angles; and the
    Euler factors in ``U3(t, p, l) = Rz(p) Ry(t) Rz(l)`` reverse as
    ``U3(-t, -l, -p)``.  The reference compiler is applied only after this
    exact logical inverse has been constructed.
    """
    gates, q0s, q1s, angles, phis, lams, width = require_kernel_arrays(kernel)
    out_gates, out_q0, out_q1, out_angles, out_phis, out_lams = [list(reversed(values)) for values in (gates, q0s, q1s, angles, phis, lams)]
    for i, gate in enumerate(out_gates):
        if gate in ROTATION_GATES:
            out_angles[i] = -out_angles[i]
        elif gate == U3:
            # U3(t, p, l) = Rz(p) Ry(t) Rz(l), so the inverse is U3(-t, -l, -p).
            out_angles[i], out_phis[i], out_lams[i] = -out_angles[i], -out_lams[i], -out_phis[i]
    return build_kernel(out_gates, out_q0, out_q1, out_angles, width, out_phis, out_lams)


@cudaq.kernel
def _mirror_kernel(width: int, kernel_a: Callable[[cudaq.qview], None], kernel_b: Callable[[cudaq.qview], None], kernel_c: Callable[[cudaq.qview], None], kernel_d: Callable[[cudaq.qview], None]):
    """Run an M1 or M2 definite-outcome mirror in temporal order.

    The paper writes products in operator order, so Eq. (34)
    ``L' g_ref(C^dag) g(C) L`` executes here as ``L`` then ``g(C)`` then
    ``g_ref(C^dag)`` then ``L'``.  Replacing ``g(C)`` by an independently
    randomized ``g_ref(C)`` gives M2, Eq. (35).  Ideally the middle pair is the
    identity and ``L' L = P``, so measurement has a known target bit string.
    """
    q = cudaq.qvector(width)
    kernel_a(q)
    kernel_b(q)
    kernel_c(q)
    kernel_d(q)
    mz(q)


BR_kernel = _mirror_kernel
RR_kernel = _mirror_kernel


@cudaq.kernel
def REF_kernel(width: int, kernel_a: Callable[[cudaq.qview], None], kernel_b: Callable[[cudaq.qview], None]):
    """Run the M3 cap-and-SPAM reference ``L' L`` of Eq. (36).

    M3 contains no QUOPS circuit.  Dividing by its mean effective polarization
    removes, to MCFE's approximation, the degradation caused by the random
    Clifford caps and native state preparation and measurement.  M3 depends on
    width but not size, so hardware experiments may reuse it across all shapes
    of the same width.
    """
    q = cudaq.qvector(width)
    kernel_a(q)
    kernel_b(q)
    mz(q)


# ---------------------------------------------------------------------------
# QUOPS circuit ensemble
# ---------------------------------------------------------------------------
#
# Supplementary Section II B defines an ensemble, not one fixed circuit.  For
# every sampled benchmark layer, the matching of active qubits, CNOT directions,
# Pauli rotation axes, and continuous angles are fresh random variables.  This
# randomness probes a broad set of computations and prevents a device from
# scoring well by specializing to a small deterministic circuit family.
#
# Public ``depth`` counts the paper's even number of elementary layers.
# The private helpers below retain benchmark-layer-pair counts for legacy
# archives; the notebook API converts paper depth to pairs exactly once.
# Result tables expose both conventions to distinguish them unambiguously.


def final_layer_qubits(width: int, zeta: float) -> int:
    """Return the number of active qubits in a legal fractional final layer.

    Supplementary Section II B restricts ``zeta`` to
    ``{1 / w, 2 / w, ..., 1}``.  Earlier versions rounded arbitrary values of
    ``zeta * w``; that silently changed the requested circuit distribution.

    The final filling fraction gives fine size resolution between complete
    benchmark layers.  It selects exactly ``zeta * w`` active qubits, then the
    same random-pairing and random-rotation rule is applied to that subset.
    """
    width = _integer_value(width, "width")
    if width < 1:
        raise ValueError("width must be at least one")
    zeta = float(zeta)
    if not math.isfinite(zeta):
        raise ValueError("zeta must be finite")

    scaled = zeta * width
    active = int(round(scaled))
    tolerance = 1.0e-12 * max(1, width)
    if not math.isclose(scaled, active, rel_tol=0.0, abs_tol=tolerance):
        raise ValueError(f"zeta={zeta} is not legal for width {width}; zeta * width must be an integer")
    if not 1 <= active <= width:
        raise ValueError(f"zeta={zeta} gives {active} final-layer qubits; must be in 1..{width}")
    return active


def paper_depth_from_benchmark_depth(depth: int) -> int:
    """Convert a legacy benchmark-layer count to the paper's depth.

    One benchmark layer is a pair of elementary layers (CNOT then rotations),
    so Supplementary Section II B and Eqs. (9)--(10) use ``d = 2 * depth``.
    """
    depth = _integer_value(depth, "benchmark-layer depth")
    if depth < 1:
        raise ValueError("benchmark-layer depth must be at least one")
    return 2 * depth


def benchmark_depth_from_paper_depth(depth: int) -> int:
    """Convert the paper's positive even elementary depth to benchmark layers."""
    depth = _integer_value(depth, "paper depth")
    if depth < 2 or depth % 2:
        raise ValueError("paper depth must be a positive even integer")
    return depth // 2


def _shape_size_from_layer_pairs(width: int, depth: int, zeta: float) -> int:
    """Return the QUOPS size ``s`` for a legacy ``(w, depth, zeta)`` input.

    ``depth`` counts benchmark layers (CNOT/rotation pairs), whereas the
    paper's ``d`` counts elementary layers and equals ``2 * depth``.  This is
    Eqs. (9)--(10) of Supplementary Section II B after that substitution::

        s_bulk  = (2 * floor(w / 2) + w) * (depth - 1)
        s_final = 2 * floor(zeta * w / 2) + zeta * w

    The factors of two count each CNOT as two operations, not as two physical
    gates.  For odd active-qubit counts, one uniformly random qubit is idle in
    the CNOT layer but still receives its independently sampled rotation.
    """
    width = _integer_value(width, "width")
    depth = _integer_value(depth, "benchmark-layer depth")
    paper_depth_from_benchmark_depth(depth)  # validates the legacy depth
    final = final_layer_qubits(width, zeta)
    bulk = (2 * (width // 2) + width) * (depth - 1)
    return bulk + 2 * (final // 2) + final


def _sample_quops_arrays_from_layer_pairs(width: int, depth: int, zeta: float, seed=None):
    """Sample a QUOPS circuit as flat gate arrays.

    Returns ``(gates, q0, q1, angles, size)``.  The legacy ``depth`` argument
    counts benchmark-layer pairs, so the paper's elementary depth is
    ``2 * depth``.  Each benchmark layer applies
    CNOTs over a uniformly random pairing of the active qubits (leaving out one
    uniformly random qubit when the count is odd) followed by ``R_P(theta)`` on
    every active qubit, with ``P`` uniform on ``{X, Y, Z}`` and ``theta``
    uniform on ``[0, 2 pi)``.  Only the final layer is partially filled.

    A random permutation of the active set is split into ordered consecutive
    pairs.  This samples an undirected matching uniformly and also makes the
    control/target orientation of every pair uniform.  If the active count is
    odd, the unpaired last element is uniformly distributed over the set.
    """
    width = _integer_value(width, "width")
    depth = _integer_value(depth, "benchmark-layer depth")
    paper_depth_from_benchmark_depth(depth)
    final_active = final_layer_qubits(width, zeta)

    # One generator supplies all choices, but every call below draws fresh
    # pseudorandom variates: active subset, pairing/orientation, axes, angles.
    rng = np.random.default_rng(seed)
    gates, q0, q1, angles = [], [], [], []
    for layer in range(depth):
        # Bulk layers act on all w qubits.  Only the last benchmark layer uses
        # the filling fraction zeta, exactly as specified in Section II B.
        if layer < depth - 1:
            active = list(range(width))
        else:
            active = sorted(int(qubit) for qubit in rng.choice(width, size=final_active, replace=False))
        # Consecutive entries of a uniform permutation form the random CNOT
        # matching.  Any final unpaired entry is the uniformly omitted qubit.
        perm = list(rng.permutation(active))
        for i in range(0, len(perm) - 1, 2):
            gates.append(CX)
            q0.append(int(perm[i]))
            q1.append(int(perm[i + 1]))
            angles.append(0.0)
        # Each active qubit receives an i.i.d. R_P(theta): P is uniform over
        # X/Y/Z and theta is uniform on [0, 2 pi), including neither endpoint
        # twice because numpy's upper bound is exclusive.
        for qubit in active:
            gates.append(ROTATION_GATES[int(rng.integers(0, 3))])
            q0.append(int(qubit))
            q1.append(0)
            angles.append(float(rng.uniform(0.0, 2.0 * math.pi)))

    return gates, q0, q1, angles, circuit_size(gates)


def _sample_quops_circuit_from_layer_pairs(width: int, depth: int, zeta: float, seed=None):
    """Sample a QUOPS circuit and return it as a built CUDA-Q kernel."""
    gates, q0, q1, angles, _ = _sample_quops_arrays_from_layer_pairs(width, depth, zeta, seed)
    return build_kernel(gates, q0, q1, angles, width)


def legal_zetas(width: int) -> tuple[float, ...]:
    """Return the admissible filling fractions ``{1/w, 2/w, ..., 1}``."""
    width = _integer_value(width, "width")
    if width < 1:
        raise ValueError("width must be at least one")
    return tuple(k / width for k in range(1, width + 1))


def _enumerate_shapes_from_layer_pairs(width_range: Sequence[int], zeta: float | None = None, *, max_depth: int = 10_000, allow_truncated: bool = False) -> list[dict]:
    """List shape dictionaries inside the utility cone.

    ``depth`` is retained as the benchmark-layer count for compatibility.
    Every dictionary also includes the unambiguous aliases
    ``benchmark_depth`` and ``paper_depth`` (the even elementary-layer depth).

    ``width_range`` is ``(min_width, max_width)``, inclusive at both ends.

    The cone is ``w**2 <= s <= w**3`` (Supplementary Eq. (19)).  With ``zeta``
    left as ``None`` every legal filling fraction is swept, which is the point
    of ``zeta``: it fills in the sizes between the multiples of ``2w`` that
    varying depth alone can reach.  Passing a float restricts the sweep to that
    single filling fraction.  ``max_depth`` is a safety cap; by default the
    function raises before doing work if that cap would silently truncate the
    requested cone.  Set ``allow_truncated=True`` only for an explicitly
    censored exploratory schedule.

    The utility cone is a reporting convention rather than a physical law.  It
    focuses the headline score on sizes large enough to be nontrivial
    (``s >= w**2``) but not so deep that width is sacrificed merely to inflate
    operation count (``s <= w**3``).  Measurements outside the cone can still
    be scientifically useful; they simply do not contribute to Q.
    """
    if len(width_range) != 2:
        raise ValueError("width_range must contain exactly (min_width, max_width)")
    low = _integer_value(width_range[0], "minimum width")
    high = _integer_value(width_range[1], "maximum width")
    if low < 1:
        raise ValueError("widths must be at least one")
    if high < low:
        raise ValueError("maximum width must not be below minimum width")
    if not isinstance(allow_truncated, (bool, np.bool_)):
        raise ValueError("allow_truncated must be a boolean")
    allow_truncated = bool(allow_truncated)
    max_depth = _integer_value(max_depth, "max_depth")
    if max_depth < 1:
        raise ValueError("max_depth must be at least one")
    shapes = []
    for width in range(low, high + 1):
        zetas = legal_zetas(width) if zeta is None else (float(zeta),)
        bulk_size = 2 * (width // 2) + width
        required_depth = max((width**3 - _shape_size_from_layer_pairs(width, 1, value)) // bulk_size + 1 for value in zetas)
        if not allow_truncated and max_depth < required_depth:
            raise ValueError(f"max_depth={max_depth} truncates the width-{width} utility cone; use at least {required_depth} or set allow_truncated=True")
        seen: dict[int, dict] = {}
        for depth in range(1, max_depth + 1):
            sizes_this_depth = []
            for value in zetas:
                size = _shape_size_from_layer_pairs(width, depth, value)
                sizes_this_depth.append(size)
                if width**2 <= size <= width**3 and size not in seen:
                    seen[size] = {"width": width, "depth": depth, "benchmark_depth": depth, "paper_depth": paper_depth_from_benchmark_depth(depth), "zeta": value, "size": size}
            if sizes_this_depth and min(sizes_this_depth) > width**3:
                break
        shapes.extend(seen[size] for size in sorted(seen))
    return shapes


# ===========================================================================
# Randomized compilation: the MCFE reference compiler g_ref
# ===========================================================================
#
# MCFE requires a compiler whose effective error is approximately stochastic
# when averaged over its randomization (Supplementary Eqs. (32)--(33)).  In
# channel notation, the average noisy reference implementation must factor as
# an approximately stochastic error channel composed with the desired ideal
# unitary, for both C and C^dag.  MCFE further assumes those two reference-error
# channels have approximately equal process polarizations.
# Supplementary Section IV B 1 specifies standard randomized compilation:
# insert an independent uniformly random Pauli layer before each two-qubit gate
# layer, insert the compensating layer after it, and compose both into the
# adjacent single-qubit layers.
#
# Averaging uniformly over the Pauli frames projects a gate-independent error
# channel onto its Pauli-diagonal part: coherent Pauli-transfer-matrix
# off-diagonal terms cancel.  That algebra is exact for an ideal Clifford layer,
# but the conclusion that a real compiled circuit has a stochastic effective
# error is still a physical assumption.  The paper gives sufficient regimes:
# single-qubit errors are small relative to two-qubit errors, approximately
# gate-independent, or already predominantly stochastic.
#
# The implementation uses the paper's Pauli-frame construction for the
# represented all-to-all-CNOT circuit, folding random Paulis into adjacent
# one-qubit gates without adding bulk gates.  Only O(w) boundary Paulis can
# remain explicit when there is no adjacent one-qubit layer into which they can
# be folded.  Each call samples an independent compilation, as required for the
# M1 and M2 ensemble averages.  On restricted-connectivity or gate-dependent
# hardware, the paper instead requires an architecture-aware reference compiler
# that performs exact routing before twirling native two-qubit Clifford layers;
# this helper alone does not provide that guarantee.


def _emit_pauli(index, qubit, gates, q0s, q1s, thetas, phis, lams):
    gate = PAULI_GATE_ID[int(index)]
    if gate is not None:
        gates.append(gate)
        q0s.append(int(qubit))
        q1s.append(0)
        thetas.append(0.0)
        phis.append(0.0)
        lams.append(0.0)


def _pauli_frame_compilation(kernel, rng):
    gates, q0s, q1s, angles, phis, lams, width = require_kernel_arrays(kernel)
    n = len(gates)

    # ``frame[q]`` is the Pauli currently carried immediately to the right of
    # the next ideal operation on qubit q.  It is bookkeeping, not automatically
    # an emitted gate.  A leading Pauli is only needed on qubits whose first
    # gate is a CNOT; for
    # any other qubit the frame is set for free by its first single-qubit gate.
    # Symmetrically, the frame need not be randomized at a single-qubit gate
    # that no later CNOT on that qubit can see, which removes the trailing
    # flush for circuits that end in a rotation layer.
    first_gate_is_cnot = [False] * width
    seen = [False] * width
    cnot_follows = [False] * n
    for index, gate in enumerate(gates):
        touched = (q0s[index], q1s[index]) if gate == CX else (q0s[index],)
        for qubit in touched:
            if not seen[qubit]:
                seen[qubit] = True
                first_gate_is_cnot[qubit] = gate == CX
    pending = [False] * width
    for index in reversed(range(n)):
        gate = gates[index]
        if gate == CX:
            pending[q0s[index]] = True
            pending[q1s[index]] = True
        else:
            cnot_follows[index] = pending[q0s[index]]

    out_gates, out_q0, out_q1, out_theta, out_phi, out_lam = [], [], [], [], [], []
    frame = [PAULI_I] * width

    for qubit in range(width):
        if first_gate_is_cnot[qubit]:
            # This is the uniformly random twirl immediately before the first
            # two-qubit layer.  I is represented by no physical instruction.
            frame[qubit] = int(rng.integers(0, 4))
            _emit_pauli(frame[qubit], qubit, out_gates, out_q0, out_q1, out_theta, out_phi, out_lam)

    for index, gate in enumerate(gates):
        if gate == CX:
            control, target = q0s[index], q1s[index]
            out_gates.append(CX)
            out_q0.append(control)
            out_q1.append(target)
            out_theta.append(0.0)
            out_phi.append(0.0)
            out_lam.append(0.0)
            # Propagate the incoming frame through the Clifford exactly.  The
            # resulting frame is the compensating Pauli on the far side because
            # P_out CX P_in = CX up to an irrelevant global phase.
            frame[control], frame[target] = conjugate_pauli_through_cnot(frame[control], frame[target])
            continue

        qubit = q0s[index]
        # A rotation layer separates successive CNOT layers in every QUOPS
        # circuit.  Sample the next independent twirl only when another CNOT
        # follows this qubit; after its last CNOT, choose I so the frame closes.
        new_frame = int(rng.integers(0, 4)) if cnot_follows[index] else PAULI_I
        # Temporal matrix order is right-to-left: absorb the incoming frame on
        # the right and the newly chosen outgoing frame on the left.  Replacing
        # the three operations by one U3 preserves the ideal circuit exactly and
        # avoids assigning the randomization an artificial gate-count overhead.
        matrix = pauli_matrix(new_frame) @ single_qubit_matrix(gate, angles[index], phis[index], lams[index]) @ pauli_matrix(frame[qubit])
        theta, phi, lam = u3_angles_from_matrix(matrix)
        out_gates.append(U3)
        out_q0.append(qubit)
        out_q1.append(0)
        out_theta.append(theta)
        out_phi.append(phi)
        out_lam.append(lam)
        frame[qubit] = new_frame

    for qubit in range(width):
        # A residual frame exists only at a boundary with no one-qubit gate to
        # absorb it.  Flushing it makes the compiled unitary exactly equal to C.
        _emit_pauli(frame[qubit], qubit, out_gates, out_q0, out_q1, out_theta, out_phi, out_lam)

    return build_kernel(out_gates, out_q0, out_q1, out_theta, width, out_phi, out_lam)


def randomized_compilation(kernel, seed=None, style: str = "pauli_frame"):
    """Pauli-frame compile a represented all-to-all-CNOT circuit.

    This is a valid MCFE ``g_ref`` when the represented CNOTs are the native
    two-qubit Clifford layers and the effective errors satisfy the stochastic
    assumptions in Supplementary Eqs. (32)--(33).  Backends that require
    routing or have gate-dependent coherent errors need an architecture-aware
    reference compiler or a separate MCFE-validity justification.

    Reusing an integer seed reproduces one compilation. Passing an advancing
    NumPy Generator draws fresh random values on each call, as required when
    averaging over the compiler expectation in Eqs. (32)--(33).
    """
    if style != "pauli_frame":
        raise ValueError(f"style must be 'pauli_frame', got {style!r}")
    return _pauli_frame_compilation(kernel, np.random.default_rng(seed))


# ===========================================================================
# MCFE caps: the random local-Clifford layers that bracket a mirror circuit
# ===========================================================================
#
# Supplementary Section IV A 1 defines three mirror circuit ensembles:
#
#     M1 = L' g_ref(c^dag) g(c)     L      (Eq. 34)
#     M2 = L' g_ref(c^dag) g_ref(c) L      (Eq. 35)
#     M3 = L'                       L      (Eq. 36)
#
# ``L`` is a layer of independent uniformly random single-qubit Cliffords, and
# ``L'`` is distributed so that ``L' L`` is an independent uniformly random
# Pauli on each qubit.  That makes every mirror circuit a definite-outcome
# circuit whose target bit string is ``P |0...0>``.
#
# The caps do two jobs.  First, the local Clifford average supplies the
# randomization used by MCFE to relate a mirror's effective polarization to a
# process polarization.  Second, choosing a random P makes the ideal output
# known without forcing it always to be all-zero.  X and Y flip |0>, whereas I
# and Z do not, which determines the target bit on each measured qubit.
#
# Rather than appending ``P`` as extra physical gates, ``L' = P L^dag`` is
# compiled into a single canonical one-qubit Clifford per qubit, so both caps
# cost exactly one gate per qubit.


# Euler triples (theta, phi, lambda) for the 24 one-qubit Cliffords under
# U3(theta, phi, lambda) = Rz(phi) Ry(theta) Rz(lambda).
_HPI = math.pi / 2.0
SINGLE_QUBIT_CLIFFORD_EULER = ((0.0, 0.0, 0.0), (_HPI, 0.0, _HPI), (_HPI, _HPI, -math.pi), (math.pi, -math.pi, 0.0), (_HPI, -math.pi, -_HPI), (_HPI, -_HPI, math.pi), (math.pi, -_HPI, -_HPI), (_HPI, -math.pi, _HPI), (_HPI, -_HPI, 0.0), (0.0, 0.0, math.pi), (_HPI, 0.0, -_HPI), (_HPI, _HPI, 0.0), (_HPI, 0.0, math.pi), (_HPI, _HPI, -_HPI), (0.0, 0.0, _HPI), (_HPI, math.pi, -math.pi), (_HPI, -_HPI, _HPI), (math.pi, math.pi, -_HPI), (_HPI, -math.pi, 0.0), (_HPI, -_HPI, -_HPI), (math.pi, 3.0 * _HPI, -math.pi), (_HPI, 0.0, 0.0), (_HPI, _HPI, -3.0 * _HPI), (0.0, 0.0, -_HPI))

NUM_ONE_QUBIT_CLIFFORDS = len(SINGLE_QUBIT_CLIFFORD_EULER)

_CLIFFORD_MATRICES = tuple(u3_matrix(*euler) for euler in SINGLE_QUBIT_CLIFFORD_EULER)


def _matching_clifford_index(unitary: np.ndarray) -> int:
    overlaps = [abs(np.trace(candidate.conj().T @ unitary)) / 2.0 for candidate in _CLIFFORD_MATRICES]
    index = int(np.argmax(overlaps))
    if overlaps[index] < 1.0 - 1e-10:
        raise RuntimeError("unitary does not match a one-qubit Clifford")
    return index


# For every sampled L index and Pauli P, compile L' = P L^dag as one canonical
# Clifford.  This preserves L' L = P up to global phase without adding P as a
# separate physical gate to the final cap.
_FINAL_CAP_CLIFFORD_INDEX = tuple(tuple(_matching_clifford_index(pauli_matrix(pauli) @ clifford.conj().T) for pauli in range(4)) for clifford in _CLIFFORD_MATRICES)


def build_clifford_layer(clifford_indices: Iterable[int]):
    """Build a tensor product of indexed one-qubit Cliffords, one gate per qubit.

    Each Clifford is emitted as a single ``u3``, including the identity, so the
    cap costs exactly one gate location per qubit and the M3 reference circuit
    measures the SPAM error of a genuine one-qubit-gate layer.

    The 24 single-qubit Clifford elements are represented canonically only up
    to global phase, which cannot affect any measurement probability or process
    channel.  Emitting the identity Clifford as U3(0, 0, 0) keeps physical cap
    depth independent of which Clifford was sampled.
    """
    indices = tuple(_integer_value(index, "one-qubit Clifford index") for index in clifford_indices)
    if not indices:
        raise ValueError("a Clifford layer must contain at least one qubit")
    if any(not 0 <= index < NUM_ONE_QUBIT_CLIFFORDS for index in indices):
        raise ValueError("one-qubit Clifford indices must be in 0..23")

    # Read the Euler-angle columns in qubit order, retaining one U3 per qubit.
    width = len(indices)
    thetas, phis, lams = zip(*(SINGLE_QUBIT_CLIFFORD_EULER[index] for index in indices))
    layer = build_kernel([U3] * width, range(width), [0] * width, thetas, width, phis, lams)
    layer.clifford_indices = indices
    return layer


def target_bitstring_from_paulis(pauli_indices: Iterable[int]) -> str:
    """Return the q0-first string produced when ``P`` acts on ``|0...0>``.

    CUDA-Q reports these sampled kernels in q0-first order, which the runtime
    tests verify explicitly.  Hamming distances in Eq. (39) are invariant to a
    common reversal, but matching each observed bit to the corresponding target
    bit is essential when the target is not symmetric.
    """
    indices = tuple(_integer_value(index, "Pauli index") for index in pauli_indices)
    if not indices:
        raise ValueError("a Pauli layer must contain at least one qubit")
    if any(not 0 <= index < 4 for index in indices):
        raise ValueError("Pauli indices must use 0=I, 1=Z, 2=X, 3=Y")
    return "".join("1" if index in (2, 3) else "0" for index in indices)


@dataclass(frozen=True)
class MCFECaps:
    """Compiled MCFE caps and the definite outcome they induce.

    ``initial`` is L and ``final`` is L' = P L^dag.  Therefore the noiseless
    cap-only circuit applies L' L = P up to global phase and deterministically
    measures ``target``.  The stored indices make this identity auditable.
    """

    initial: Callable[[cudaq.qview], None]
    final: Callable[[cudaq.qview], None]
    target: str
    clifford_indices: tuple[int, ...]
    final_clifford_indices: tuple[int, ...]
    pauli_indices: tuple[int, ...]

    @property
    def width(self) -> int:
        return len(self.clifford_indices)


def build_mcfe_caps(clifford_indices: Iterable[int], pauli_indices: Iterable[int]) -> MCFECaps:
    """Build ``L``, ``L' = P L^dag``, and the resulting target bit string."""
    clifford_indices = tuple(_integer_value(index, "one-qubit Clifford index") for index in clifford_indices)
    pauli_indices = tuple(_integer_value(index, "Pauli index") for index in pauli_indices)
    if len(clifford_indices) != len(pauli_indices):
        raise ValueError("Clifford and Pauli layers must have the same width")
    target = target_bitstring_from_paulis(pauli_indices)
    initial = build_clifford_layer(clifford_indices)
    final_clifford_indices = tuple(_FINAL_CAP_CLIFFORD_INDEX[clifford][pauli] for clifford, pauli in zip(clifford_indices, pauli_indices))
    final = build_clifford_layer(final_clifford_indices)
    return MCFECaps(initial=initial, final=final, target=target, clifford_indices=clifford_indices, final_clifford_indices=final_clifford_indices, pauli_indices=pauli_indices)


def sample_mcfe_caps(width: int, seed=None) -> MCFECaps:
    """Sample independent uniform local Cliffords and target Paulis for MCFE.

    Independence is required across qubits and mirror instances.  Reusing one
    cap for M1, M2, and M3 would create correlations not represented by the
    estimator or its cluster bootstrap, so the notebook samples three
    separate caps on every mirror iteration.
    """
    width = _integer_value(width, "MCFE cap width")
    if width < 1:
        raise ValueError("MCFE cap width must be at least one")
    rng = np.random.default_rng(seed)
    cliffords = tuple(int(value) for value in rng.integers(0, NUM_ONE_QUBIT_CLIFFORDS, size=width))
    paulis = tuple(int(value) for value in rng.integers(0, 4, size=width))
    return build_mcfe_caps(cliffords, paulis)


# ===========================================================================
# Noise models
# ===========================================================================
#
# This helper supplies a controlled simulator model for examples and regression
# tests; it is not part of the QUOPS protocol.  A hardware benchmark measures
# the backend's native noisy implementation and normally passes no synthetic
# noise model.  Applying the same synthetic channel to every emitted gate lets
# the tests compare MCFE with independently calculated Choi-state process
# polarizations under a known CPTP model.


# Every single-qubit gate the QUOPS kernels can emit.  ``u3`` matters because
# the pauli_frame reference compiler emits general single-qubit gates; omitting
# it would leave the reference circuits partly noiseless and bias MCFE upwards.
ONE_QUBIT_GATE_KEYS = ("rx", "ry", "rz", "x", "y", "z", "u3")
TWO_QUBIT_GATE_KEYS = ("cx",)


def make_depolarizing_noise(p_1q: float = 0.0, p_2q: float = 0.0):
    """Build a uniform CUDA-Q depolarizing model over every emitted gate.

    CUDA-Q's probability here is the total probability of applying a nonidentity
    Pauli.  Consequently a 1Q channel has process polarization ``1 - 4p/3`` and
    a 2Q channel has ``1 - 16p/15``; ``p`` itself is not the process
    polarization parameter used by the paper's alternative depolarizing-channel
    convention.  Covering U3 and boundary X/Y/Z gates is essential because the
    reference compiler and MCFE caps emit them even when C contains only axial
    rotations and CNOTs.
    """
    p_1q = _noise_probability("p_1q", p_1q)
    p_2q = _noise_probability("p_2q", p_2q)
    noise = cudaq.NoiseModel()
    if p_1q:
        channel = cudaq.DepolarizationChannel(float(p_1q))
        for gate in ONE_QUBIT_GATE_KEYS:
            noise.add_all_qubit_channel(gate, channel)
    if p_2q:
        for gate in TWO_QUBIT_GATE_KEYS:
            noise.add_all_qubit_channel(gate, cudaq.Depolarization2(float(p_2q)))
    return noise


def make_helios1_noise(gate_error_1q: float = 2.50e-5, gate_error_2q: float = 7.90e-4, spam_0: float = 8.00e-4, spam_1: float = 1.60e-4):
    """Approximate the supplied Helios-1 snapshot dated 2025-11-05.

    Gate errors are benchmark infidelities, converted to nonidentity-Pauli
    probabilities for a leakage-free depolarizing surrogate. RX/RY/X/Y/H/U3
    receive one 1Q channel; virtual Z rotations receive none. One 2Q channel
    is assigned per abstract CX, without native decomposition overhead.

    ``spam_0`` and ``spam_1`` combine preparation and measurement error. They
    are represented once, as P(1|0) and P(0|1) at terminal Z measurement;
    this is not a model of mid-circuit measurement and reset. Crosstalk,
    leakage, transport and memory noise are not separately reconstructed.

    Source: user-supplied Quantinuum Nexus Helios-1 calibration snapshot,
    calibration date 2025-11-05 (PDF export 2026-09-10). Reported values and
    uncertainties: 1Q infidelity 2.50e-5 +/- 1.00e-6; 2Q infidelity 7.90e-4
    +/- 2.00e-5; SPAM 0 8.00e-4 +/- 1.00e-4; SPAM 1 1.60e-4 +/- 5.00e-5.
    Measurement crosstalk was 4.80e-5 +/- 1.00e-6; memory error was unspecified.
    Calibration uncertainties are provenance, not extra channel probabilities.

    For dimension d and average infidelity r, the nonidentity-Pauli probability
    is p=(d+1)*r/d: the defaults give p_1q=3.75e-5 and p_2q=9.875e-4. U3 is
    treated as one physical rotation with virtual Z phases. Native Helios CX
    decomposition also needs two 1Q rotations; their added cost is omitted.
    The model is a reproducible surrogate, not the native emulator or a hardware
    performance prediction. Definitions and native-gate mapping are documented at:
    https://docs.quantinuum.com/systems/user_guide/hardware_user_guide/performance_validation.html
    https://docs.quantinuum.com/systems/user_guide/hardware_user_guide/helios.html#native-gate-set
    Use a noise-capable target such as ``density-matrix-cpu``.
    """
    gate_error_1q = _noise_probability("gate_error_1q", gate_error_1q)
    gate_error_2q = _noise_probability("gate_error_2q", gate_error_2q)
    spam_0 = _noise_probability("spam_0", spam_0)
    spam_1 = _noise_probability("spam_1", spam_1)
    # For dimension d, average infidelity r = d*p/(d+1), not p itself.
    p_1q = _noise_probability("1.5 * gate_error_1q", 1.5 * gate_error_1q)
    p_2q = _noise_probability("1.25 * gate_error_2q", 1.25 * gate_error_2q)
    noise = cudaq.NoiseModel()
    if p_1q:
        for gate in ("rx", "ry", "x", "y", "h", "u3"):
            noise.add_all_qubit_channel(gate, cudaq.DepolarizationChannel(p_1q))
    if p_2q:
        noise.add_all_qubit_channel("cx", cudaq.Depolarization2(p_2q))
    if spam_0 or spam_1:
        # CUDA-Q applies this channel before mz: K0 preserves, K1 flips 0, K2 flips 1.
        preserve = np.diag(np.sqrt([1.0 - spam_0, 1.0 - spam_1])).astype(complex)
        flip_zero = np.array([[0.0, 0.0], [np.sqrt(spam_0), 0.0]], dtype=complex)
        flip_one = np.array([[0.0, np.sqrt(spam_1)], [0.0, 0.0]], dtype=complex)
        noise.add_all_qubit_channel("mz", cudaq.KrausChannel([preserve, flip_zero, flip_one]))
    return noise


# ===========================================================================
# Statistics: estimators, confidence bounds, and scores
#
# Process fidelity compares the noisy channel with its ideal unitary (SI
# Eq. (5)); process polarization rescales it (Eq. (6)). The benchmark target is
# the ensemble mean process polarization, not mean output success probability.
# Effective polarization is the Hamming-distance statistic of a definite-outcome
# mirror (Eq. (39)); MCFE combines its ensemble means to infer the target.
# The standard threshold is gamma_0 = 1/sqrt(e), per SI Eqs. (16)--(17).
# Normal and percentile confidence procedures below remain approximate under
# the MCFE assumptions. The complete predeclared scan uses Bonferroni shares.
# ===========================================================================

QUOPS_THRESHOLD = 1.0 / math.sqrt(math.e)
DEFAULT_FAMILYWISE_ALPHA = 0.05
MIN_BOOTSTRAP_TAIL_DRAWS = 10
MCFE_ESTIMATORS = ("ratio_of_means", "mean_of_ratios")


def _validate_estimator(estimator: str) -> str:
    estimator = str(estimator)
    if estimator not in MCFE_ESTIMATORS:
        raise ValueError(f"estimator must be one of {MCFE_ESTIMATORS}, got {estimator!r}")
    return estimator


def hamming_distance(a: str, b: str) -> int:
    """Return the Hamming distance between equal-length bit strings."""
    a = str(a).replace(" ", "")
    b = str(b).replace(" ", "")
    if len(a) != len(b):
        raise ValueError(f"bit strings must have the same length: {a!r}, {b!r}")
    return sum(x != y for x, y in zip(a, b))


def effective_polarization_from_counts(counts, width: int, target: str | None = None) -> float:
    """Estimate a mirror's effective polarization using Eq. (39).

    Let ``h_k`` be the observed frequency of outcomes at Hamming distance ``k``
    from the mirror's known target.  First compute

        moment = sum_k (-1/2)**k h_k,

    then affinely rescale it as

        lambda_hat = (moment - 4**(-w)) / (1 - 4**(-w)).

    An ideal definite-outcome circuit has ``lambda_hat = 1``.  Uniform random
    output has expected moment ``4**(-w)`` and therefore polarization zero.
    Unlike ordinary success probability, errors at every Hamming distance
    contribute with alternating, exponentially decreasing weight.  This is the
    statistic for which the local-Clifford MCFE averaging theory applies.
    """
    width = _positive_integer("width", width)
    target = "0" * width if target is None else str(target).replace(" ", "")
    if len(target) != width or set(target) - {"0", "1"}:
        raise ValueError(f"target must be a {width}-bit binary string")
    count_items = counts.items() if hasattr(counts, "items") else counts.to_dict().items() if hasattr(counts, "to_dict") else dict(counts).items()
    normalized_counts = {}
    for bitstring, value in count_items:
        bitstring = str(bitstring).replace(" ", "")
        if len(bitstring) != width or set(bitstring) - {"0", "1"}:
            raise ValueError(f"count key {bitstring!r} must be a {width}-bit binary string")
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("shot counts must be non-negative integers") from exc
        if not math.isfinite(numeric) or numeric < 0.0 or not numeric.is_integer():
            raise ValueError("shot counts must be non-negative integers")
        normalized_counts[bitstring] = normalized_counts.get(bitstring, 0) + int(numeric)
    total = sum(normalized_counts.values())
    if total <= 0:
        raise ValueError("counts must contain at least one shot")
    # Grouping explicitly by k is unnecessary: weighting each observed outcome
    # by (-1/2)^Hamming-distance is algebraically the same empirical sum.
    moment = sum(count * (-0.5) ** hamming_distance(bitstring, target) for bitstring, count in normalized_counts.items()) / total
    # This inverse-dimension form is Eq. (39) divided through by 4**w.  It is
    # mathematically identical but remains finite for hundreds of qubits.
    inverse_dimension_sq = 4.0 ** (-width)
    return float((moment - inverse_dimension_sq) / (1.0 - inverse_dimension_sq))


def normalized_mcfe_polarization(gamma_br: float, gamma_rr: float, gamma_ref: float, *, min_denominator: float = 1.0e-12) -> float:
    """Combine mean M1, M2, and M3 effective polarizations using Eq. (40).

    In the notebook's names, ``br`` is M1, ``rr`` is M2, and ``ref`` is M3:

        gamma(C) approximately lambda_M1 / sqrt(lambda_M2 * lambda_M3).

    The square-root normalization removes the reference-compiler contribution
    measured by M2 and the cap/SPAM contribution measured by M3.  The relation
    is an MCFE approximation under Eqs. (32)--(33), not an exact identity for an
    arbitrary noise process.  A nonpositive or numerically tiny denominator is
    physically unstable, so it returns NaN and downstream inference fails closed.
    """
    gamma_br = float(gamma_br)
    gamma_rr = float(gamma_rr)
    gamma_ref = float(gamma_ref)
    min_denominator = float(min_denominator)
    if not math.isfinite(min_denominator) or min_denominator < 0.0:
        raise ValueError("min_denominator must be finite and non-negative")
    if not all(np.isfinite(value) for value in (gamma_br, gamma_rr, gamma_ref)) or gamma_rr <= min_denominator or gamma_ref <= min_denominator:
        return float("nan")
    return float(gamma_br / (math.sqrt(gamma_rr) * math.sqrt(gamma_ref)))


def _strict_mean(values) -> float:
    array = np.asarray(tuple(values), dtype=float)
    if array.size == 0 or not np.isfinite(array).all():
        return float("nan")
    return float(array.mean())


def mcfe_polarization(br_by_circuit, rr_by_circuit, ref_samples, estimator: str = "ratio_of_means") -> float:
    """Estimate shape polarization using Eq. (53) or circuit-wise Eqs. (56)--(57).

    ``ratio_of_means`` is the paper's lower-variance Eq. (53) estimator, which
    omits the covariance term in Eqs. (58)--(59) and can therefore be slightly
    optimistic.  ``mean_of_ratios`` retains the circuit-wise normalization but
    generally has higher variance.  Both retain MCFE's physical assumptions.

    Data layout is ``[circuit][mirror]`` for M1 and M2, and a separate flat M3
    sample.  With equal mirror counts, Eq. (53) first averages all M1 values and
    all M2 values, then forms one ratio.  Eqs. (56)--(57) instead normalize each
    circuit's M1 mean by its own M2 mean before averaging over circuits.  It is
    scientifically invalid to pair arbitrary individual M1, M2, and M3 rows and
    average those row-wise ratios; those rows are independent ensemble draws.

    Equation (54) permits M2 to be omitted only in the special case where the
    implementation under test is itself the same randomized reference compiler,
    ``g == g_ref``.  This general estimator deliberately retains M2 because the
    notebook benchmarks the unrandomized implementation ``g(C)`` against a
    separately randomized ``g_ref(C)``.
    """
    estimator = _validate_estimator(estimator)
    if len(br_by_circuit) != len(rr_by_circuit) or not br_by_circuit or any(not group for group in (*br_by_circuit, *rr_by_circuit)) or len({len(group) for group in br_by_circuit}) != 1 or len({len(group) for group in rr_by_circuit}) != 1:
        return float("nan")
    gamma_ref = _strict_mean(ref_samples)
    if not np.isfinite(gamma_ref) or gamma_ref <= 0.0:
        return float("nan")
    if estimator == "ratio_of_means":
        # Equation (53): lower sampling variance, but it drops the circuit-level
        # covariance term written explicitly in Eqs. (58)--(59).
        gamma_br = _strict_mean(value for group in br_by_circuit for value in group)
        gamma_rr = _strict_mean(value for group in rr_by_circuit for value in group)
        return normalized_mcfe_polarization(gamma_br, gamma_rr, gamma_ref)
    # Equations (56)--(57): preserve the M1/M2 normalization within each shared
    # target circuit C_k, then average the resulting circuit estimates.
    br_means = tuple(_strict_mean(values) for values in br_by_circuit)
    rr_means = tuple(_strict_mean(values) for values in rr_by_circuit)
    if not all(np.isfinite(value) for value in (*br_means, *rr_means)):
        return float("nan")
    return _strict_mean(normalized_mcfe_polarization(gamma_br, gamma_rr, gamma_ref) for gamma_br, gamma_rr in zip(br_means, rr_means))


@dataclass(frozen=True)
class ShapePolarizationData:
    """Effective-polarization observations for one shape ``(w, s)``.

    ``br[k][l]`` and ``rr[k][l]`` are M1 and M2 observations conditional on the
    kth independently sampled QUOPS circuit.  M3 is independent of C and is
    stored once as ``ref[l]``; it may contain more observations than one
    circuit's mirror count.  Preserving this hierarchy is necessary both for
    Eq. (53) and for statistically coherent cluster resampling.

    ``depth`` is the paper's positive even elementary-layer count ``d``.
    ``paper_depth`` is an alias; ``benchmark_depth`` counts the ``d // 2``
    CNOT/rotation pairs. These units match all public circuit and scan helpers.

    ``metadata`` may hold operational quantities such as M1/M2 elapsed time and
    total M1/M2 shots.  Those values do not enter the polarization estimator;
    they enter the separate rate calculation.
    """
    width: int
    size: int
    depth: int
    zeta: float
    br: tuple[tuple[float, ...], ...]
    rr: tuple[tuple[float, ...], ...]
    ref: tuple[float, ...]
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        width = _positive_integer("width", self.width)
        size = _positive_integer("size", self.size)
        depth = _integer_value(self.depth, "paper depth")
        zeta = float(self.zeta)
        br = tuple(tuple(float(value) for value in group) for group in self.br)
        rr = tuple(tuple(float(value) for value in group) for group in self.rr)
        ref = tuple(float(value) for value in self.ref)
        expected_size = shape_size(width, depth, zeta)
        if size != expected_size:
            raise ValueError(f"size={size} does not match generated size {expected_size}")
        if len(br) != len(rr) or not br:
            raise ValueError("br and rr must contain the same nonzero number of circuits")
        if any(not group for group in (*br, *rr)) or not ref:
            raise ValueError("every M1/M2 circuit group and the M3 ensemble must be nonempty")
        if len({len(group) for group in br}) != 1 or len({len(group) for group in rr}) != 1:
            raise ValueError("each of M1 and M2 must use a fixed mirror count across circuits")
        if not all(np.isfinite(value) for group in (*br, *rr) for value in group):
            raise ValueError("M1/M2 polarizations must be finite")
        if not all(np.isfinite(value) for value in ref):
            raise ValueError("M3 polarizations must be finite")
        inverse_dimension_sq = 4.0 ** (-width)
        support_lower = (-0.5 - inverse_dimension_sq) / (1.0 - inverse_dimension_sq)
        if any(value < support_lower - 1.0e-12 or value > 1.0 + 1.0e-12 for group in (*br, *rr, ref) for value in group):
            raise ValueError(f"effective mirror polarizations must lie in [{support_lower}, 1]")
        object.__setattr__(self, "width", width)
        object.__setattr__(self, "size", size)
        object.__setattr__(self, "depth", depth)
        object.__setattr__(self, "zeta", zeta)
        object.__setattr__(self, "br", br)
        object.__setattr__(self, "rr", rr)
        object.__setattr__(self, "ref", ref)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def num_circuits(self) -> int:
        return len(self.br)

    @property
    def benchmark_depth(self) -> int:
        return self.depth // 2

    @property
    def paper_depth(self) -> int:
        return self.depth

    def estimate(self, estimator: str = "ratio_of_means") -> float:
        return mcfe_polarization(self.br, self.rr, self.ref, estimator)


def bootstrap_polarization(data: ShapePolarizationData, *, estimator: str = "ratio_of_means", resamples: int = 2000, seed=None) -> np.ndarray:
    """Resample whole independent circuit clusters jointly for M1 and M2.

    Each observed circuit's M1 and M2 means already contain the variability of
    its independently sampled mirrors and shots. Resampling those mirrors a
    second time would double count within-circuit variation. Joint circuit
    indices preserve the M1/M2 covariance induced by the common target circuit.
    M3 uses an independently resampled width-only reference ensemble.

    One mirror per circuit is sufficient when multiple circuits are observed.
    The resulting uncertainty is approximate and concerns the chosen MCFE
    estimator; it does not account for reference-compiler bias or model error.
    """
    estimator = _validate_estimator(estimator)
    resamples = _positive_integer("resamples", resamples)
    rng = np.random.default_rng(seed)
    draws = np.empty(resamples, dtype=float)
    br = np.asarray([np.mean(group) for group in data.br], dtype=float)
    rr = np.asarray([np.mean(group) for group in data.rr], dtype=float)
    ref = np.asarray(data.ref, dtype=float)
    for draw_index in range(resamples):
        circuit_choice = rng.integers(0, data.num_circuits, size=data.num_circuits)
        reference_mean = float(ref[rng.integers(0, ref.size, size=ref.size)].mean())
        if estimator == "ratio_of_means":
            draws[draw_index] = normalized_mcfe_polarization(float(br[circuit_choice].mean()), float(rr[circuit_choice].mean()), reference_mean)
        elif reference_mean <= 1.0e-12 or np.any(rr[circuit_choice] <= 1.0e-12):
            draws[draw_index] = float("nan")
        else:
            draws[draw_index] = float(np.mean(br[circuit_choice] / np.sqrt(rr[circuit_choice] * reference_mean)))
    return draws


def bootstrap_lower_bound(draws, *, alpha: float = DEFAULT_FAMILYWISE_ALPHA) -> float:
    """Return the empirical one-sided percentile endpoint without interpolation.

    Nonfinite ratios are retained as negative infinity. Discarding unstable
    denominators would condition on favorable draws and bias the endpoint up.
    This order statistic alone has no finite-sample coverage guarantee.
    """
    alpha = _probability("alpha", alpha)
    values = np.asarray(tuple(draws), dtype=float)
    if values.size == 0:
        return float("nan")
    values = np.where(np.isfinite(values), values, -np.inf)
    index = int(math.floor((values.size - 1) * alpha))
    return float(np.partition(values, index)[index])


def _validate_interval_method(interval_method: str) -> str:
    if interval_method not in ("normal", "percentile"):
        raise ValueError("interval_method must be 'normal' or 'percentile'")
    return interval_method


def _validate_threshold(threshold: float) -> float:
    if isinstance(threshold, (bool, np.bool_)):
        raise ValueError("threshold must be numeric, not boolean")
    threshold = float(threshold)
    if not math.isfinite(threshold) or not 0.0 < threshold <= 1.0:
        raise ValueError("threshold must be finite and in (0, 1]")
    return threshold


def validate_inference_config(family_size: int, *, familywise_alpha: float = DEFAULT_FAMILYWISE_ALPHA, bootstrap_resamples: int = 2000, interval_method: str = "normal", num_circuits: int | None = None, num_mirrors: int | None = None, num_m3_mirrors: int | None = None) -> dict:
    """Validate a predeclared Bonferroni scan before expensive acquisition.

    The normal approximation estimates a standard error, so it does not need
    empirical alpha-tail resolution. A percentile interval needs at least ten
    requested tail draws; that is a numerical guard, not a coverage theorem.
    Multiple independent circuits and M3 observations remain essential.
    """
    family_size = _positive_integer("family_size", family_size)
    familywise_alpha = _probability("familywise_alpha", familywise_alpha)
    bootstrap_resamples = _positive_integer("bootstrap_resamples", bootstrap_resamples)
    interval_method = _validate_interval_method(interval_method)
    alpha = familywise_alpha / family_size
    if bootstrap_resamples < 2:
        raise ValueError("at least two bootstrap resamples are required to estimate a standard error")
    if interval_method == "percentile" and bootstrap_resamples * alpha < MIN_BOOTSTRAP_TAIL_DRAWS:
        raise ValueError(f"bootstrap tail is under-resolved: need at least {math.ceil(MIN_BOOTSTRAP_TAIL_DRAWS / alpha)} resamples for this declared family")
    if num_circuits is not None and _positive_integer("num_circuits", num_circuits) < 2:
        raise ValueError("inference requires at least two independent circuits")
    if num_mirrors is not None:
        _positive_integer("num_mirrors", num_mirrors)
    if num_m3_mirrors is not None and _positive_integer("num_m3_mirrors", num_m3_mirrors) < 2:
        raise ValueError("inference requires at least two independent M3 observations")
    return {"family_size": family_size, "familywise_alpha": familywise_alpha, "alpha_per_shape": alpha, "bootstrap_resamples": bootstrap_resamples, "interval_method": interval_method, "confidence_is_approximate": True}


def _bounded_degenerate_lower_bound(data: ShapePolarizationData, alpha: float, estimator: str) -> float:
    """Guard a degenerate bootstrap using simultaneous Hoeffding mean bounds.

    A single Hamming statistic lies in [(-1/2-4**(-w))/(1-4**(-w)), 1].
    The same support holds for averages of mirrors within a circuit. Applying
    three one-sided bounds with alpha/3 bounds the ratio-of-population-means
    estimand when its numerator and denominators are positive. This conservative
    fallback assumes independent circuits and independent M3 observations; it
    still cannot bound the physical MCFE approximation error. It is not a bound
    for the different mean-of-ratios estimand, which remains unestablished here.
    """
    if estimator != "ratio_of_means":
        return float("-inf")
    inverse_dimension_sq = 4.0 ** (-data.width)
    support_lower = (-0.5 - inverse_dimension_sq) / (1.0 - inverse_dimension_sq)
    observations = [value for group in (*data.br, *data.rr) for value in group] + list(data.ref)
    if any(value < support_lower - 1.0e-12 or value > 1.0 + 1.0e-12 for value in observations):
        return float("-inf")
    circuit_radius = (1.0 - support_lower) * math.sqrt(math.log(3.0 / alpha) / (2.0 * data.num_circuits))
    reference_radius = (1.0 - support_lower) * math.sqrt(math.log(3.0 / alpha) / (2.0 * len(data.ref)))
    br_lower = float(np.mean(data.br)) - circuit_radius
    rr_upper = min(1.0, float(np.mean(data.rr)) + circuit_radius)
    ref_upper = min(1.0, float(np.mean(data.ref)) + reference_radius)
    return br_lower / math.sqrt(rr_upper * ref_upper) if br_lower > 0.0 and rr_upper > 0.0 and ref_upper > 0.0 else float("-inf")


def pass_statistics(data: ShapePolarizationData, *, estimator: str = "ratio_of_means", threshold: float = QUOPS_THRESHOLD, alpha: float = DEFAULT_FAMILYWISE_ALPHA, resamples: int = 2000, seed=None, interval_method: str = "normal") -> dict:
    """Test a shape with an explicitly approximate MCFE lower confidence bound.

    The default follows the current paper's Methods normal approximation,
    centered on the observed estimate: estimate - z_(1-alpha) * bootstrap_SE.
    The standard error comes from whole-circuit cluster resampling. Percentile
    endpoints remain available explicitly. Neither construction has general
    finite-sample coverage, and Bonferroni cannot repair a miscalibrated test.

    A zero empirical variance is not proof of a deterministic population. Such
    samples use a conservative bounded-mean fallback for ratio_of_means, or
    remain unestablished for mean_of_ratios. Denominator failures make a normal
    interval unestablished rather than being discarded. The legacy column
    bootstrap_lower is retained as an alias for polarization_lower.
    """
    from statistics import NormalDist
    estimator = _validate_estimator(estimator)
    interval_method = _validate_interval_method(interval_method)
    alpha = _probability("alpha", alpha)
    threshold = _validate_threshold(threshold)
    resamples = _positive_integer("resamples", resamples)
    point = data.estimate(estimator)
    enough_independent_data = data.num_circuits >= 2 and len(data.ref) >= 2
    enough_tail_resolution = interval_method != "percentile" or resamples * alpha >= MIN_BOOTSTRAP_TAIL_DRAWS
    draws = bootstrap_polarization(data, estimator=estimator, resamples=resamples, seed=seed) if enough_independent_data and enough_tail_resolution and resamples >= 2 and np.isfinite(point) else np.array([], dtype=float)
    finite_draws = draws[np.isfinite(draws)]
    bootstrap_stderr = float(finite_draws.std(ddof=1)) if finite_draws.size > 1 else float("nan")
    lower = float("nan")
    fallback = "none"
    if draws.size:
        if finite_draws.size != draws.size and interval_method == "normal":
            lower = float("-inf")
        elif finite_draws.size == draws.size and bootstrap_stderr <= 1.0e-14 * max(1.0, abs(point)):
            fallback = "hoeffding_component_means" if estimator == "ratio_of_means" else "degenerate_mean_of_ratios_unestablished"
            lower = _bounded_degenerate_lower_bound(data, alpha, estimator)
        elif interval_method == "normal":
            lower = point + NormalDist().inv_cdf(alpha) * bootstrap_stderr
        else:
            lower = bootstrap_lower_bound(draws, alpha=alpha)
    if not enough_independent_data:
        decision, reason = "invalid", "inference requires at least two independent circuits and two independent M3 observations"
    elif resamples < 2:
        decision, reason = "invalid", "at least two bootstrap resamples are required"
    elif not enough_tail_resolution:
        decision, reason = "invalid", f"bootstrap tail is under-resolved: resamples * alpha must be at least {MIN_BOOTSTRAP_TAIL_DRAWS}"
    elif not np.isfinite(point):
        decision, reason = "invalid", "MCFE point estimate is undefined because a normalization denominator is unstable"
    elif fallback != "none":
        decision = "approximate_pass" if lower >= threshold else "not_certified"
        reason = "zero empirical variance; conservative bounded-mean fallback clears the threshold" if decision == "approximate_pass" else "zero empirical variance cannot establish the threshold with these sample counts"
    elif not np.isfinite(lower):
        decision, reason = "not_certified", "bootstrap denominator failures prevent a reliable lower bound"
    elif lower >= threshold:
        decision, reason = "approximate_pass", f"approximate one-sided {interval_method} lower bound clears the threshold"
    else:
        decision, reason = "not_certified", f"approximate one-sided {interval_method} lower bound does not clear the threshold"
    valid_fraction = float(finite_draws.size / draws.size) if draws.size else 0.0
    return {"width": data.width, "size": data.size, "depth": data.depth, "depth_unit": "paper", "benchmark_depth": data.benchmark_depth, "paper_depth": data.paper_depth, "zeta": data.zeta, "num_circuits": data.num_circuits, "estimator": estimator, "mean_polarization": point, "threshold": threshold, "alpha": alpha, "bootstrap_lower": lower, "polarization_lower": lower, "bootstrap_stderr": bootstrap_stderr, "bootstrap_resamples": resamples, "bootstrap_valid": int(finite_draws.size), "bootstrap_valid_fraction": valid_fraction, "confidence_method": f"cluster_bootstrap_{interval_method}", "interval_method": interval_method, "interval_fallback": fallback, "confidence_is_approximate": True, "passes_threshold": bool(np.isfinite(point) and point >= threshold), "decision": decision, "decision_reason": reason, "passes": decision == "approximate_pass"}


# ===========================================================================
# QUOPS scores and operational rates
# ===========================================================================

def utility_cone_bounds(width: int) -> tuple[int, int]:
    """Return the inclusive size bounds ``w**2 <= s <= w**3`` of Eq. (19).

    The cone determines eligibility for the headline Q score, not whether a
    circuit can be generated or measured.  Results outside it may be reported as
    diagnostics but must not enlarge Q.
    """
    width = _positive_integer("width", width)
    return width**2, width**3


def in_utility_cone(width: int, size: int) -> bool:
    """Return whether a shape lies in the utility cone."""
    low, high = utility_cone_bounds(width)
    return low <= _integer_value(size, "size") <= high


def quops_score_from_summary(summary_df: pd.DataFrame, *, pass_column: str = "passes") -> tuple[int, tuple[int, int] | None]:
    """Return Eq. (20)'s largest passing tested size and its ``(w, s)`` shape.

    This function only summarizes rows actually supplied.  A finite scan that
    ends while its largest shape still passes gives a right-censored lower bound
    on capability, not evidence that the true boundary occurs at that size.
    Statistical multiplicity control must already be encoded by how those rows
    were scheduled and tested.
    """
    # In the paper's gated analysis, region-only passes cannot increase the score.
    if pass_column == "passes" and "score_passes" in summary_df:
        pass_column = "score_passes"
    required = {"width", "size", pass_column}
    if summary_df.empty or not required.issubset(summary_df.columns):
        return 0, None
    eligible = summary_df[pass_column].fillna(False).eq(True)
    if "decision" in summary_df:
        eligible &= ~summary_df["decision"].eq("invalid")
    if "measured" in summary_df:
        eligible &= summary_df["measured"].fillna(False).eq(True)
    if "mean_polarization" in summary_df:
        eligible &= np.isfinite(pd.to_numeric(summary_df["mean_polarization"], errors="coerce"))
    passed = summary_df[eligible]
    if passed.empty:
        return 0, None
    passed = passed[[in_utility_cone(row.width, row.size) for row in passed[["width", "size"]].itertuples(index=False)]]
    if passed.empty:
        return 0, None
    best = passed.sort_values(["size", "width"], ascending=[False, True]).iloc[0]
    best_width = _positive_integer("best width", best["width"])
    best_size = _positive_integer("best size", best["size"])
    return best_size, (best_width, best_size)


def quops_rate(size: int, gamma: float, tau_wall: float, m1_m2_shots: int, *, kept_shots: int | None = None) -> float:
    """Return Eq. (72), including the optional postselection shot factor.

    The formula is

        Omega = 2 s gamma_bar**2 sqrt(N_total * N_kept) / tau_wall.

    ``s * gamma_bar**2`` is the paper's useful-operation weighting.  The factor
    two compensates for an M1/M2 mirror being approximately twice the base
    circuit depth.  ``m1_m2_shots`` is the total number of completed executions
    of all M1 and M2 mirror circuits used for the shape. ``kept_shots`` defaults
    to that total for backwards compatibility. When postselection is supplied,
    gamma must be estimated on accepted data under a justified protocol; this
    rate helper does not define an acceptance rule or validate that protocol.

    For an operational hardware result, ``tau_wall`` should cover acquisition
    of all M1/M2 data, including unavoidable execution-system overhead, while
    excluding cloud queue and data-transfer time.  Width-only M3 acquisition may
    be excluded because the same reference data can be reused across sizes.
    """
    size = _nonnegative_integer("size", size)
    m1_m2_shots = _nonnegative_integer("m1_m2_shots", m1_m2_shots)
    kept_shots = m1_m2_shots if kept_shots is None else _nonnegative_integer("kept_shots", kept_shots)
    if kept_shots > m1_m2_shots:
        raise ValueError("kept_shots cannot exceed total m1_m2_shots")
    if isinstance(tau_wall, (bool, np.bool_)):
        raise ValueError("tau_wall must be finite and positive")
    tau_wall = float(tau_wall)
    if not math.isfinite(tau_wall) or tau_wall <= 0.0:
        raise ValueError("tau_wall must be finite and positive")
    if isinstance(gamma, (bool, np.bool_)):
        raise ValueError("gamma used for a physical rate must be numeric, not boolean")
    gamma = float(gamma)
    if not math.isfinite(gamma):
        return float("nan")
    if not 0.0 <= gamma <= 1.0:
        raise ValueError("gamma used for a physical rate must lie in [0, 1]")
    return float(2.0 * size * gamma**2 * math.sqrt(m1_m2_shots * kept_shots) / tau_wall)


def omega_simulated(size: int, gamma: float, tau_wall: float, m1_m2_shots: int, *, kept_shots: int | None = None) -> float:
    """Evaluate Eq. (72) with a local simulator timing proxy.

    The minimal notebook sums the elapsed spans of individual M1/M2 CUDA-Q
    calls.  That includes work inside those calls but excludes M3, circuit
    construction, analysis, idle gaps, and some system overhead.  Host load,
    compilation caches, and JIT warm-up can all affect it, so this value is a
    reproducibility diagnostic and must not be presented as hardware Omega.
    """
    return quops_rate(size, gamma, tau_wall, m1_m2_shots, kept_shots=kept_shots)


# ===========================================================================
# Scan scheduling and shape-level summaries
# ===========================================================================

def _default_shape_schedule_from_layer_pairs(shapes: Iterable[dict] | Sequence[int]) -> list[dict]:
    """Validate and order the complete finite family before acquisition.

    The paper permits a user-specified shape-selection algorithm and a strategy
    for allocating the 5% familywise error budget.  This deterministic order is
    used by the minimal notebook to acquire the complete predeclared family.
    ``summarize_quops_runs`` applies a Bonferroni allocation, so the result does
    not depend on adaptive stopping. Per-shape coverage remains approximate.

    Duplicate ``(w,s)`` hypotheses are rejected to prevent accidental retesting
    and incorrect allocation of the declared familywise error budget.
    """
    shapes = list(shapes)
    if shapes and not isinstance(shapes[0], dict):
        shapes = _enumerate_shapes_from_layer_pairs(shapes)
    canonical = []
    seen = set()
    for candidate in shapes:
        if not isinstance(candidate, dict):
            raise ValueError("every scheduled shape must be a dictionary")
        if "width" not in candidate or "depth" not in candidate:
            raise ValueError("every scheduled shape requires width and depth")
        if candidate.get("depth_unit", "layer_pairs") != "layer_pairs":
            raise ValueError("legacy archive schedules require depth in layer_pairs; use default_shape_schedule for paper depth")
        width = _positive_integer("shape width", candidate["width"])
        depth = _positive_integer("shape depth", candidate["depth"])
        zeta = float(candidate.get("zeta", 1.0))
        size = _shape_size_from_layer_pairs(width, depth, zeta)
        if "size" in candidate and _integer_value(candidate["size"], "scheduled size") != size:
            raise ValueError(f"scheduled size {candidate['size']} does not match computed size {size}")
        hypothesis = (width, size)
        if hypothesis in seen:
            raise ValueError(f"duplicate scheduled hypothesis (width={width}, size={size})")
        seen.add(hypothesis)
        canonical.append({**candidate, "width": width, "depth": depth, "benchmark_depth": depth, "paper_depth": paper_depth_from_benchmark_depth(depth), "zeta": zeta, "size": size})
    ordered = sorted(canonical, key=lambda shape: (shape["size"], shape["width"]))
    for schedule_index, shape in enumerate(ordered):
        shape["schedule_index"] = schedule_index
    return ordered


def _summarize_quops_runs_from_layer_pairs(run_df: pd.DataFrame, *, familywise_alpha: float = DEFAULT_FAMILYWISE_ALPHA, bootstrap_resamples: int = 2000, seed=None, estimator: str = "ratio_of_means", threshold: float = QUOPS_THRESHOLD, interval_method: str = "normal", declared_shapes: Iterable[dict] | None = None, family_size: int | None = None, adjust_familywise: bool = True) -> pd.DataFrame:
    """Aggregate mirror records using Eq. (53) and a declared Bonferroni family.

    M1/M2 are clustered by the independently drawn circuit. The row-wise gamma
    diagnostic is never averaged to estimate a shape. Both shape estimators and
    their difference are reported to expose sensitivity to circuit covariance;
    that difference is a diagnostic, not a correction or bound on MCFE bias.

    Supply declared_shapes from the acquisition schedule to preserve the family
    after filtering or missing acquisition. Missing declared shapes are invalid
    rows. An explicit family_size also preserves the alpha budget, but cannot
    identify missing shapes. Omitting both retains the legacy observed-row
    behavior, which is labeled observed_rows_only and is unsuitable for
    reconstructing a predeclared family from a filtered CSV.

    Optional m1_m2_kept_shots enables Eq. (72); omitted means all shots were kept.
    Timings from local simulator calls remain throughput proxies. Empty input
    is accepted with an explicit declaration so missing data cannot become a
    smaller testing family. Integer fields are validated before conversion.
    """
    familywise_alpha = _probability("familywise_alpha", familywise_alpha)
    estimator = _validate_estimator(estimator)
    threshold = _validate_threshold(threshold)
    interval_method = _validate_interval_method(interval_method)
    bootstrap_resamples = _positive_integer("bootstrap_resamples", bootstrap_resamples)
    explicit_family_size = family_size is not None
    declared = _default_shape_schedule_from_layer_pairs(declared_shapes) if declared_shapes is not None else None
    if declared is not None and not declared:
        raise ValueError("declared_shapes must contain at least one shape")
    if family_size is not None:
        family_size = _positive_integer("family_size", family_size)
    if declared is not None and family_size is not None and family_size != len(declared):
        raise ValueError("family_size must equal the complete declared_shapes count")
    required = {"shape", "width", "depth", "size", "zeta", "circuit_index", "mirror_index", "gamma_br", "gamma_rr", "gamma_ref", "br_rr_seconds", "m1_m2_shots"}
    if run_df.empty and (declared is not None or family_size is not None):
        data = run_df.reindex(columns=sorted(required | set(run_df.columns))).copy()
    else:
        if not required.issubset(run_df.columns):
            raise ValueError(f"run_df is missing required columns: {sorted(required - set(run_df.columns))}")
        data = run_df.copy()
    positive_integer_columns = ("width", "depth", "size", "m1_m2_shots")
    nonnegative_integer_columns = ("circuit_index", "mirror_index") + (("m1_m2_kept_shots",) if "m1_m2_kept_shots" in data.columns else ())
    for column in positive_integer_columns:
        data[column] = pd.Series([_positive_integer(column, value) for value in data[column]], index=data.index, dtype="int64")
    for column in nonnegative_integer_columns:
        data[column] = pd.Series([_nonnegative_integer(column, value) for value in data[column]], index=data.index, dtype="int64")
    numeric_columns = ("zeta", "gamma_br", "gamma_rr", "gamma_ref", "br_rr_seconds")
    for column in numeric_columns:
        if any(isinstance(value, (bool, np.bool_)) for value in data[column]):
            raise ValueError(f"{column} must be numeric, not boolean")
        data[column] = pd.to_numeric(data[column], errors="coerce")
    if data[list(numeric_columns)].isna().any().any() or not np.isfinite(data[list(numeric_columns)].to_numpy(dtype=float)).all():
        raise ValueError("run_df contains missing or non-finite required measurements")
    if (data["br_rr_seconds"] <= 0.0).any():
        raise ValueError("every run needs positive BR/RR time")
    if "m1_m2_kept_shots" in data.columns and (data["m1_m2_kept_shots"] > data["m1_m2_shots"]).any():
        raise ValueError("m1_m2_kept_shots cannot exceed total m1_m2_shots")
    if data.duplicated(["width", "size", "circuit_index", "mirror_index"]).any():
        raise ValueError("run_df contains duplicate shape/circuit/mirror rows")
    grouped = list(data.groupby(["width", "depth", "size", "zeta"], sort=True))
    observed_shapes = {(int(width), int(size)): (int(depth), float(zeta)) for (width, depth, size, zeta), _group in grouped}
    if len(observed_shapes) != len(grouped):
        raise ValueError("each (width, size) hypothesis must have a single depth and zeta")
    if declared is not None:
        planned_shapes = {(shape["width"], shape["size"]): (shape["depth"], shape["zeta"]) for shape in declared}
        if any(key not in planned_shapes or planned_shapes[key][0] != value[0] or not math.isclose(planned_shapes[key][1], value[1], rel_tol=1.0e-12, abs_tol=1.0e-12) for key, value in observed_shapes.items()):
            raise ValueError("observed shapes do not match declared_shapes")
        family_size = len(declared)
    family_size = len(grouped) if family_size is None else family_size
    if family_size < len(grouped):
        raise ValueError("family_size cannot be smaller than the observed hypothesis count")
    config = validate_inference_config(family_size, familywise_alpha=familywise_alpha, bootstrap_resamples=bootstrap_resamples, interval_method=interval_method)
    shape_alpha = config["alpha_per_shape"] if adjust_familywise else familywise_alpha
    config["alpha_per_shape"] = shape_alpha
    schedule_indices = {(shape["width"], shape["size"]): shape["schedule_index"] for shape in declared or ()}
    provenance = "declared_shapes" if declared is not None else "explicit_family_size" if explicit_family_size else "observed_rows_only"
    summaries = []
    for shape_index, ((width, depth, size, zeta), group) in enumerate(grouped):
        ordered_groups = [circuit_group.sort_values("mirror_index") for _, circuit_group in group.groupby("circuit_index", sort=True)]
        br = tuple(tuple(circuit_group["gamma_br"].astype(float)) for circuit_group in ordered_groups)
        rr = tuple(tuple(circuit_group["gamma_rr"].astype(float)) for circuit_group in ordered_groups)
        ref = tuple(group.sort_values(["circuit_index", "mirror_index"])["gamma_ref"].astype(float))
        shape_data = ShapePolarizationData(width=width, size=size, depth=paper_depth_from_benchmark_depth(depth), zeta=zeta, br=br, rr=rr, ref=ref)
        inference_index = schedule_indices.get((width, size), shape_index)
        stats = pass_statistics(shape_data, estimator=estimator, threshold=threshold, alpha=shape_alpha, resamples=bootstrap_resamples, seed=None if seed is None else _integer_value(seed, "seed") + 1_000_003 * (inference_index + 1), interval_method=interval_method)
        stats["depth"] = depth
        stats.pop("depth_unit", None)
        stats["shape"] = str(group.iloc[0]["shape"])
        stats["num_mirrors"] = len(ordered_groups[0])
        stats["mean_gamma_br"] = float(group["gamma_br"].mean())
        stats["mean_gamma_rr"] = float(group["gamma_rr"].mean())
        stats["mean_gamma_ref"] = float(group["gamma_ref"].mean())
        stats["ratio_of_means"] = shape_data.estimate("ratio_of_means")
        stats["mean_of_ratios"] = shape_data.estimate("mean_of_ratios")
        stats["estimator_gap"] = stats["ratio_of_means"] - stats["mean_of_ratios"]
        stats["tau_wall_seconds"] = float(group["br_rr_seconds"].sum())
        stats["m1_m2_shots"] = sum(int(value) for value in group["m1_m2_shots"])
        stats["m1_m2_kept_shots"] = sum(int(value) for value in group["m1_m2_kept_shots"]) if "m1_m2_kept_shots" in group.columns else stats["m1_m2_shots"]
        stats["omega_simulated"] = omega_simulated(size, stats["mean_polarization"], stats["tau_wall_seconds"], stats["m1_m2_shots"], kept_shots=stats["m1_m2_kept_shots"]) if np.isfinite(stats["mean_polarization"]) and 0.0 <= stats["mean_polarization"] <= 1.0 else np.nan
        stats["measured"] = True
        summaries.append(stats)
    for shape in declared or ():
        if (shape["width"], shape["size"]) not in observed_shapes:
            summaries.append({**shape, "shape": f"w={shape['width']},s={shape['size']}", "num_circuits": 0, "num_mirrors": 0, "estimator": estimator, "mean_polarization": np.nan, "threshold": threshold, "alpha": shape_alpha, "bootstrap_lower": np.nan, "polarization_lower": np.nan, "bootstrap_stderr": np.nan, "bootstrap_resamples": bootstrap_resamples, "confidence_method": f"cluster_bootstrap_{interval_method}", "interval_method": interval_method, "interval_fallback": "none", "confidence_is_approximate": True, "passes_threshold": False, "decision": "invalid", "decision_reason": "declared shape has no acquired observations", "passes": False, "measured": False})
    columns = sorted({key for row in summaries for key in row} | {"width", "size", "passes", "mean_polarization", "decision"})
    result = pd.DataFrame(summaries, columns=columns)
    result["multiplicity"] = "bonferroni" if adjust_familywise else "unadjusted"
    result["familywise_alpha"] = familywise_alpha
    result["family_size"] = family_size
    result["observed_family_size"] = len(grouped)
    result["family_declaration"] = provenance
    result.attrs.update({**config, "observed_family_size": len(grouped), "family_declaration": provenance, "declared_shapes": declared})
    return result.sort_values(["size", "width"]).reset_index(drop=True)


# ===========================================================================
# Public API: the paper's even elementary-layer depth d
# ===========================================================================


def shape_size(width: int, depth: int, zeta: float) -> int:
    """Return Eqs. (9)--(10)'s QUOPS size for paper depth d (a positive even integer).

    Every pair of elementary layers contains CNOTs then single-qubit rotations.
    Thus s = (2 floor(w/2) + w) (d/2 - 1) + 2 floor(zeta*w/2) + zeta*w.
    The final active-qubit count zeta*w must be an integer in 1..w.
    """
    return _shape_size_from_layer_pairs(width, benchmark_depth_from_paper_depth(depth), zeta)


def sample_quops_arrays(width: int, depth: int, zeta: float, seed=None):
    """Sample flat gate arrays for paper depth d: d/2 CNOT/rotation layer pairs."""
    return _sample_quops_arrays_from_layer_pairs(width, benchmark_depth_from_paper_depth(depth), zeta, seed)


def sample_quops_circuit(width: int, depth: int, zeta: float, seed=None):
    """Build a CUDA-Q circuit with paper depth d before mirroring or compilation."""
    return _sample_quops_circuit_from_layer_pairs(width, benchmark_depth_from_paper_depth(depth), zeta, seed)


def _paper_shape_to_layer_pairs(shape: dict) -> dict:
    """Validate paper-unit metadata before using the legacy inference internals."""
    if not isinstance(shape, dict) or "depth" not in shape:
        raise ValueError("every paper-depth shape requires depth")
    depth = _integer_value(shape["depth"], "paper depth")
    layer_pairs = benchmark_depth_from_paper_depth(depth)
    if shape.get("depth_unit", "paper") != "paper":
        raise ValueError("public helpers require depth_unit='paper'; convert legacy layer-pair data explicitly")
    for name, expected in (("paper_depth", depth), ("benchmark_depth", layer_pairs)):
        if name in shape and _integer_value(shape[name], name) != expected:
            raise ValueError(f"{name}={shape[name]} conflicts with paper depth {depth}")
    return {**shape, "depth": layer_pairs, "benchmark_depth": layer_pairs, "paper_depth": depth, "depth_unit": "layer_pairs"}


def _paper_shape_from_layer_pairs(shape: dict) -> dict:
    layer_pairs = _integer_value(shape["depth"], "benchmark-layer depth")
    depth = paper_depth_from_benchmark_depth(layer_pairs)
    return {**shape, "depth": depth, "benchmark_depth": layer_pairs, "paper_depth": depth, "depth_unit": "paper"}


def enumerate_shapes(width_range: Sequence[int], zeta: float | None = None, *, max_depth: int = 20_000, allow_truncated: bool = False) -> list[dict]:
    """Enumerate utility-cone shapes; depth and max_depth use the paper's even d.

    Width bounds are inclusive. With zeta=None, enumerate all legal final-layer
    filling fractions. An insufficient max_depth raises unless truncation is
    explicitly allowed. Each shape also records its benchmark_depth=d/2.
    """
    shapes = _enumerate_shapes_from_layer_pairs(width_range, zeta, max_depth=benchmark_depth_from_paper_depth(max_depth), allow_truncated=allow_truncated)
    return [_paper_shape_from_layer_pairs(shape) for shape in shapes]


def default_shape_schedule(shapes: Iterable[dict] | Sequence[int]) -> list[dict]:
    """Validate the predeclared family using the paper's positive even depth d.

    A two-element width range instead enumerates the utility cone. Explicit
    benchmark_depth/paper_depth aliases must agree with depth; legacy marked
    schedules are rejected so saved layer-pair depths cannot be reinterpreted.
    """
    candidates = list(shapes)
    if candidates and not isinstance(candidates[0], dict):
        candidates = enumerate_shapes(candidates)
    scheduled = _default_shape_schedule_from_layer_pairs([_paper_shape_to_layer_pairs(shape) for shape in candidates])
    return [_paper_shape_from_layer_pairs(shape) for shape in scheduled]


def summarize_quops_runs(run_df: pd.DataFrame, *, familywise_alpha: float = DEFAULT_FAMILYWISE_ALPHA, bootstrap_resamples: int = 2000, seed=None, estimator: str = "ratio_of_means", threshold: float = QUOPS_THRESHOLD, interval_method: str = "normal", declared_shapes: Iterable[dict] | None = None, family_size: int | None = None, score_sequence: Sequence[tuple[int, int]] | None = None) -> pd.DataFrame:
    """Summarize measurements whose depth is the paper's even d.

    Polarization, confidence bounds and rate use the existing MCFE inference.
    Only depth metadata is converted internally. Returned rows and declared
    shapes in attrs retain paper units. Save depth_unit='paper' with raw rows;
    old archives use private layer-pair entry points to retain their meaning.

    Supply score_sequence as predeclared (width, size) pairs to use the paper's
    gated score and Hochberg capability tests (SI IV D 1--3). The default 0.05
    applies separately to both families. Omitting it retains archive behavior.
    """
    if run_df.attrs.get("depth_unit", "paper") != "paper":
        raise ValueError("measurements require depth_unit='paper'; convert legacy layer-pair data explicitly")
    data = run_df.copy()
    if "depth" in data:
        normalized = [_paper_shape_to_layer_pairs(row) for row in data.to_dict("records")]
        for name in ("depth", "benchmark_depth", "paper_depth"):
            data[name] = pd.Series([row[name] for row in normalized], index=data.index, dtype="int64")
        data["depth_unit"] = "layer_pairs"
    data.attrs["depth_unit"] = "layer_pairs"
    declared = default_shape_schedule(declared_shapes) if declared_shapes is not None else None
    if score_sequence is not None:
        if declared is None or interval_method != "normal" or estimator != "ratio_of_means":
            raise ValueError("the paper workflow requires declared_shapes, normal intervals and the ratio-of-means estimator")
        familywise_alpha = _probability("familywise_alpha", familywise_alpha)
        if familywise_alpha >= 0.5:
            raise ValueError("each of the paper's two family error budgets must be below 0.5")
        score_sequence = _validate_score_sequence(score_sequence, declared)
    normalized_declaration = [_paper_shape_to_layer_pairs(shape) for shape in declared] if declared is not None else None
    result = _summarize_quops_runs_from_layer_pairs(data, familywise_alpha=familywise_alpha, bootstrap_resamples=bootstrap_resamples, seed=seed, estimator=estimator, threshold=threshold, interval_method=interval_method, declared_shapes=normalized_declaration, family_size=family_size, adjust_familywise=score_sequence is None)
    if "depth" in result:
        result["benchmark_depth"] = result["depth"].astype("int64")
        result["depth"] = result["benchmark_depth"].map(paper_depth_from_benchmark_depth)
        result["paper_depth"] = result["depth"]
    result["depth_unit"] = "paper"
    result.attrs.update({"depth_unit": "paper", "declared_shapes": declared})
    return _apply_paper_quops_tests(result, score_sequence, familywise_alpha) if score_sequence is not None else result


def _validate_score_sequence(sequence, declared_shapes):
    """Check SI Eqs. (63)--(64); the caller must choose S before inspecting data."""
    sequence = [(_positive_integer("score width", width), _positive_integer("score size", size)) for width, size in sequence]
    declared = {(shape["width"], shape["size"]) for shape in declared_shapes}
    if not sequence or any(shape not in declared or not in_utility_cone(*shape) for shape in sequence):
        raise ValueError("score_sequence must contain declared shapes inside the utility cone")
    if any(right[1] <= left[1] for left, right in zip(sequence, sequence[1:])):
        raise ValueError("score_sequence sizes must be strictly increasing in their predeclared order")
    return sequence


def _apply_paper_quops_tests(summary, score_sequence, alpha):
    """Apply the paper's normal test, gated score, and Hochberg region tests.

    SI Eqs. (60)--(68): S stops at its first failed test; R excludes every
    member of S, including its untested tail. Missing/undefined tests cannot
    pass or reduce the declared family. A zero standard error is undefined
    in Eq. (61), so it is marked invalid instead of assigned a certain pass.
    The normal and MCFE approximations, independence/positive dependence for
    Hochberg, and monotonicity for downward closure remain assumptions.
    """
    from statistics import NormalDist
    result = summary.copy()
    sequence = _validate_score_sequence(score_sequence, result[["width", "size"]].to_dict("records"))
    indices = {(int(row["width"]), int(row["size"])): index for index, row in result.iterrows()}
    valid = result["measured"].fillna(False).eq(True) & result["decision"].ne("invalid") & np.isfinite(result["mean_polarization"]) & np.isfinite(result["bootstrap_stderr"]) & result["bootstrap_stderr"].gt(0.0) & result.get("bootstrap_valid_fraction", pd.Series(np.nan, index=result.index)).eq(1.0)
    result["p_value"] = np.nan
    result["interval_fallback"] = "none"
    for index in result.index[valid]:
        row = result.loc[index]
        result.loc[index, "p_value"] = NormalDist().cdf((row["threshold"] - row["mean_polarization"]) / row["bootstrap_stderr"])
        result.loc[index, ["polarization_lower", "bootstrap_lower"]] = row["mean_polarization"] - NormalDist().inv_cdf(1.0 - alpha) * row["bootstrap_stderr"]
    result["adjusted_p_value"] = np.nan
    result["test_family"] = "Hochberg (R)"
    result["tested"] = False
    result["score_passes"] = False
    result["passes"] = False
    result["decision"] = "not_tested"
    result["decision_reason"] = "after the first unsuccessful gated test"
    result.loc[~valid, ["polarization_lower", "bootstrap_lower"]] = np.nan

    # Each reached member of S uses the full 5% budget; stop on the first non-pass.
    gate_open = True
    for shape in sequence:
        index = indices[shape]
        result.loc[index, "test_family"] = "Gated score (S)"
        if not gate_open:
            result.loc[index, ["p_value", "polarization_lower", "bootstrap_lower"]] = np.nan
            continue
        passed = bool(valid.loc[index] and result.loc[index, "p_value"] <= alpha)
        result.loc[index, ["tested", "score_passes", "passes"]] = [True, passed, passed]
        result.loc[index, "decision"] = "pass" if passed else "fail" if valid.loc[index] else "invalid"
        result.loc[index, "decision_reason"] = "normal one-sided test at the gated significance level" if valid.loc[index] else "normal test undefined: missing/insufficient data, unstable denominator or zero standard error"
        gate_open = passed

    # Eq. (67): step up over all of R; invalid/missing cases keep their places with p=1.
    remaining = result.index[result["test_family"].eq("Hochberg (R)")]
    if len(remaining):
        ordered = result.loc[remaining, "p_value"].fillna(1.0).sort_values(kind="stable")
        scaled = ordered.to_numpy() * np.arange(len(ordered), 0, -1)
        adjusted = np.minimum(1.0, np.minimum.accumulate(scaled[::-1])[::-1])
        result.loc[ordered.index, "adjusted_p_value"] = adjusted
        result.loc[remaining, "tested"] = True
        for index in remaining:
            passed = bool(valid.loc[index] and result.loc[index, "adjusted_p_value"] <= alpha)
            result.loc[index, "passes"] = passed
            result.loc[index, "decision"] = "pass" if passed else "fail" if valid.loc[index] else "invalid"
            result.loc[index, "decision_reason"] = "Hochberg-adjusted p-value compared with the region-family budget" if valid.loc[index] else "normal test undefined: missing/insufficient data, unstable denominator or zero standard error"

    result["multiplicity"] = "gated_hochberg"
    result["score_confidence"] = 1.0 - alpha
    result["region_confidence"] = 1.0 - 2.0 * alpha
    result.attrs.update({"multiplicity": "gated_hochberg", "score_sequence": sequence, "score_confidence": 1.0 - alpha, "region_confidence": 1.0 - 2.0 * alpha, "score": quops_score_from_summary(result)[0], "region_assumption": "mean polarization is non-increasing in width and size"})
    return result


# ===========================================================================
# Plotting: capability, polarization, timing and rate
# ===========================================================================

def _power_of_two_ticks(min_value: float, max_value: float) -> list[float]:
    lower = int(math.floor(math.log2(max(float(min_value), 1.0))))
    upper = int(math.ceil(math.log2(max(float(max_value), 1.0))))
    return [float(2**exponent) for exponent in range(lower, upper + 1)]


def _power_of_two_label(value: float, _position=None) -> str:
    if value <= 0 or not np.isfinite(value):
        return ""
    exponent = round(math.log2(value))
    return f"$2^{{{exponent}}}$" if abs(value - 2**exponent) <= max(1.0e-9, 1.0e-9 * value) else ""


def _paper_tick_label(value, _position=None):
    """Use ordinary decade labels at small scales, as in the paper's Fig. 2a-c."""
    return f"{value:,.0f}" if 1 <= value < 100_000 else f"$10^{{{round(math.log10(value))}}}$"


def _paper_plot_axes(ax, *, numeric_ticks=False):
    """Apply the main paper's decade scales, fine ticks and uncluttered axes."""
    from matplotlib.ticker import FuncFormatter, LogFormatterMathtext, LogLocator, NullFormatter
    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10, numticks=8))
        axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
        axis.set_major_formatter(FuncFormatter(_paper_tick_label) if numeric_ticks else LogFormatterMathtext(base=10))
        axis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", direction="out", length=4, width=0.8, labelsize=10, color="#555555")
    ax.tick_params(which="minor", direction="out", length=2, width=0.6, color="#777777")
    for spine in ax.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(0.8)
    ax.grid(False, which="both")
    ax.set_axisbelow(True)


def plot_quops_scan(summary_df: pd.DataFrame, *, title: str = "QUOPS capability", label: str | None = None, pass_column: str = "passes", threshold: float | None = None, jitter_seed: int = 123, shape_figsize: tuple[float, float] = (9.0, 4.8), polarization_figsize: tuple[float, float] = (9.0, 4.2), dpi: int = 140, ax=None):
    """Plot measured capability points and polarization-versus-size diagnostics.

    Main-paper Figs. 1d and 2a-c supply the decade size/width axes, filled/open
    circles and a distinct score diamond, shown here in NVIDIA green shades.
    The pale green utility cone is geometric. With score_sequence results,
    overlay Eq. (68)'s downward closure under the paper's monotonicity assumption.
    Legacy summaries show only the tested points.
    Open circles mean the threshold was not established, not proven failure.
    Gray squares distinguish invalid or unmeasured shapes. Coordinates are
    exact; jitter_seed is retained for compatibility but no jitter is applied.

    The second figure's symmetric error bars are bootstrap standard deviations
    for visualization.  Classification uses the separate one-sided lower bound,
    so those error bars must not be read as the decision interval.
    Pass ``ax`` to draw the capability plot in an existing subplot.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    required = {"width", "size", pass_column, "mean_polarization"}
    if not required.issubset(summary_df.columns):
        raise ValueError(f"summary_df is missing required columns: {sorted(required - set(summary_df.columns))}")
    threshold_values = pd.to_numeric(summary_df["threshold"], errors="coerce").dropna().unique() if "threshold" in summary_df.columns else np.array([])
    if threshold is None and len(threshold_values) > 1:
        raise ValueError("a scan plot requires a common threshold")
    threshold = _validate_threshold(float(threshold_values[0]) if threshold is None and len(threshold_values) == 1 else QUOPS_THRESHOLD if threshold is None else threshold)
    if len(threshold_values) and not np.allclose(threshold_values, threshold, rtol=0.0, atol=1.0e-12):
        raise ValueError("plot threshold must match the threshold used for shape decisions")
    plot_df = summary_df.copy().sort_values(["width", "size"])
    for column in ("width", "size", "mean_polarization"):
        plot_df[column] = pd.to_numeric(plot_df[column], errors="coerce")
    error_column = "bootstrap_stderr" if "bootstrap_stderr" in plot_df.columns else "stderr"
    plot_df["yerr"] = pd.to_numeric(plot_df[error_column], errors="coerce").fillna(0.0) if error_column in plot_df.columns else 0.0
    plot_df = plot_df.dropna(subset=["width", "size"])
    for row in plot_df[["width", "size"]].itertuples(index=False):
        _positive_integer("width", row.width)
        _positive_integer("size", row.size)
    if plot_df.empty:
        raise ValueError("summary_df does not contain plottable rows")
    invalid_mask = plot_df["decision"].eq("invalid") if "decision" in plot_df.columns else pd.Series(False, index=plot_df.index)
    if "measured" in plot_df:
        invalid_mask |= ~plot_df["measured"].fillna(False).eq(True)
    plot_df["inference_invalid"] = invalid_mask | ~np.isfinite(plot_df["mean_polarization"])
    plot_df[pass_column] = plot_df[pass_column].fillna(False).eq(True) & ~plot_df["inference_invalid"]
    score, best_shape = quops_score_from_summary(plot_df, pass_column=pass_column)
    paper_protocol = "score_passes" in plot_df
    x_min, x_max = 1.0, 10.0 ** max(1, math.ceil(math.log10(float(plot_df["size"].max()))))
    x_padding = 10.0**0.02
    y_min, y_max = 1.0 / x_padding, 10.0 ** max(1, math.ceil(math.log10(float(plot_df["width"].max())))) * x_padding
    x_grid = np.geomspace(x_min, x_max, 500)
    cone_lower = x_grid ** (1.0 / 3.0)
    cone_upper = x_grid**0.5
    standalone = ax is None
    if standalone:
        fig_shape, ax = plt.subplots(figsize=shape_figsize, dpi=dpi)
    else:
        fig_shape = ax.figure
    _paper_plot_axes(ax, numeric_ticks=True)
    ax.fill_between(x_grid, cone_lower, cone_upper, color="#EDF5DD", linewidth=0, label=r"Utility cone: $w^2\leq s\leq w^3$", zorder=0).set_gid("quops-utility-cone")
    ax.plot(x_grid, cone_lower, color="#9DBB6B", linewidth=0.8, zorder=1)
    ax.plot(x_grid, cone_upper, color="#9DBB6B", linewidth=0.8, zorder=1)
    pass_df = plot_df[plot_df[pass_column]]
    not_tested = plot_df["decision"].eq("not_tested") if "decision" in plot_df else pd.Series(False, index=plot_df.index)
    fail_df = plot_df[~plot_df[pass_column] & ~plot_df["inference_invalid"] & ~not_tested]
    invalid_df = plot_df[plot_df["inference_invalid"]]
    if paper_protocol and not pass_df.empty:
        # Eq. (68): exact step boundary of the downward closure of declared passes.
        region_x = sorted({x_min, x_max, *pass_df["size"].astype(float)})
        region_y = [max([1.0, *pass_df.loc[pass_df["size"].ge(right), "width"].astype(float)]) for right in region_x[1:]] + [1.0]
        confidence = float(plot_df["region_confidence"].iloc[0])
        ax.fill_between(region_x, 1.0, region_y, step="post", color="#B4D984", alpha=0.4, label=f"Inferred region ({confidence:.0%}; monotonicity assumed)", zorder=0.5).set_gid("quops-capability-region")
    if not pass_df.empty:
        ax.scatter(pass_df["size"], pass_df["width"], s=34, marker="o", color="#76B900", linewidths=0.8, label="Pass" if paper_protocol else "approximate MCFE pass", zorder=3).set_gid("quops-pass")
    if not fail_df.empty:
        ax.scatter(fail_df["size"], fail_df["width"], s=34, marker="o", facecolors="white", edgecolors="#76B900", linewidths=1.1, label="threshold not established", zorder=3).set_gid("quops-unestablished")
    if not invalid_df.empty:
        ax.scatter(invalid_df["size"], invalid_df["width"], s=34, marker="s", facecolors="none", edgecolors="#777777", linewidths=1.1, label="invalid inference / unmeasured", zorder=3).set_gid("quops-invalid")
    untested_df = plot_df[not_tested & ~plot_df["inference_invalid"]]
    if not untested_df.empty:
        ax.scatter(untested_df["size"], untested_df["width"], s=34, marker="^", facecolors="none", edgecolors="#777777", linewidths=1.1, label="Not tested: gate closed", zorder=3).set_gid("quops-not-tested")
    if best_shape is not None:
        score_label = f"QUOPS score: $\\hat Q={score}$" if paper_protocol else f"Approximate tested score: $Q_{{\\mathrm{{tested}}}}={score}$"
        ax.scatter([best_shape[1]], [best_shape[0]], s=80, marker="d", facecolors="#4B7600", edgecolors="#222222", linewidths=0.8, label=score_label, zorder=4).set_gid("quops-score")
    ax.set_xlim(x_min / x_padding, x_max * x_padding)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("Circuit size, $s$ (QUOPS)", fontsize=11)
    ax.set_ylabel("Circuit width, $w$ (qubits)", fontsize=11)
    ax.set_title(title if label is None else f"{title}\n{label}", fontsize=11)
    boundary_x = x_min * (x_max / x_min)**0.22
    ax.text(boundary_x, boundary_x**0.5 * 1.10, r"$s=w^2$", fontsize=9, color="#666666")
    ax.text(boundary_x, boundary_x**(1.0 / 3.0) / 1.22, r"$s=w^3$", fontsize=9, color="#666666")
    if best_shape is None:
        ax.text(0.97, 0.04, "No gated score pass" if paper_protocol else "No approximate pass in this scan", transform=ax.transAxes, ha="right", fontsize=9, color="#555555")
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.9, fontsize=9)
    if standalone:
        fig_shape.tight_layout()
    fig_polarization, ax = plt.subplots(figsize=polarization_figsize, dpi=dpi)
    ax.set_prop_cycle(color=["#76B900", "#4B7600", "#9ACB44", "#345200", "#B4D984", "#608F20"])
    for width, group in plot_df.groupby("width"):
        group = group.sort_values("size")
        ax.errorbar(group["size"], group["mean_polarization"], yerr=group["yerr"], marker="o", markersize=4, linewidth=1.2, capsize=2, alpha=0.85, label=f"w={int(width)}")
    threshold_label = r"$1/\sqrt{e}$" if math.isclose(threshold, QUOPS_THRESHOLD, rel_tol=0.0, abs_tol=1.0e-12) else f"threshold = {threshold:g}"
    ax.axhline(threshold, color="#4B7600", linestyle="--", linewidth=1.2, label=threshold_label)
    ax.set_xscale("log", base=2)
    ax.set_xlim(x_min / x_padding, x_max * x_padding)
    ax.set_xticks(_power_of_two_ticks(x_min, x_max))
    ax.xaxis.set_major_formatter(FuncFormatter(_power_of_two_label))
    ax.set_xlabel("Size, $s$")
    ax.set_ylabel("Mean MCFE polarization")
    ax.set_title("Polarization by shape" if label is None else f"Polarization by shape\n{label}")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig_polarization.tight_layout()
    return fig_shape, fig_polarization


def plot_quops_timing(summary_df: pd.DataFrame, *, time_column: str = "tau_wall_seconds", label: str | None = None, figsize: tuple[float, float] = (9.0, 4.2), dpi: int = 140):
    """Plot the cumulative M1/M2 simulator-call time recorded for each shape.

    This is a raw simulator diagnostic with one aggregate observation per
    shape, not an operational hardware timing distribution.  It excludes M3
    reference calls and analysis, while CUDA-Q invocation overhead inside each
    timed M1/M2 call is included.  With only one aggregate timing observation
    per shape, no sampling-uncertainty error bar is implied.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    required = {"width", "size", time_column}
    if not required.issubset(summary_df.columns):
        raise ValueError(f"summary_df is missing required columns: {sorted(required - set(summary_df.columns))}")
    plot_df = summary_df.copy().sort_values(["width", "size"])
    for column in ("width", "size", time_column):
        plot_df[column] = pd.to_numeric(plot_df[column], errors="coerce")
    if plot_df.empty or plot_df[["width", "size", time_column]].isna().any().any() or not np.isfinite(plot_df[["width", "size", time_column]].to_numpy(dtype=float)).all():
        raise ValueError("summary_df contains missing or non-finite timing data")
    if (plot_df[time_column] <= 0.0).any():
        raise ValueError(f"{time_column} must be positive")
    for row in plot_df[["width", "size"]].itertuples(index=False):
        _positive_integer("width", row.width)
        _positive_integer("size", row.size)
    x_min = 2.0 ** math.floor(math.log2(max(1.0, float(plot_df["size"].min()))))
    x_max = 2.0 ** math.ceil(math.log2(max(2.0, float(plot_df["size"].max()))))
    x_max = x_min * 2.0 if x_max <= x_min else x_max
    x_padding = 2.0**0.035
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_prop_cycle(color=["#76B900", "#4B7600", "#9ACB44", "#345200", "#B4D984", "#608F20"])
    for width, group in plot_df.groupby("width"):
        group = group.sort_values("size")
        ax.plot(group["size"], group[time_column], marker="o", markersize=4, linewidth=1.2, alpha=0.85, label=f"w={int(width)}")
    ax.set_xscale("log", base=2)
    ax.set_xlim(x_min / x_padding, x_max * x_padding)
    ax.set_xticks(_power_of_two_ticks(x_min, x_max))
    ax.xaxis.set_major_formatter(FuncFormatter(_power_of_two_label))
    ax.set_xlabel("Size, $s$")
    ax.set_ylabel("Cumulative M1/M2 simulator-call time (s)")
    ax.set_title("M1/M2 simulator timing by shape" if label is None else f"M1/M2 simulator timing by shape\n{label}")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    return fig


def plot_quops_rate(summary_df: pd.DataFrame, *, rate_column: str = "omega_simulated", pass_column: str = "passes", label: str | None = None, figsize: tuple[float, float] = (9.0, 4.2), dpi: int = 140, ax=None):
    """Plot the saved shape-wise QUOPS/s values and identify the rate at Q.

    Each point is Equation (72) evaluated for one measured shape, using the
    shape-level MCFE ratio-of-means estimate rather than a mean of arbitrary
    per-mirror ratios.  The minimal notebook supplies local CUDA-Q call spans,
    so its default ``omega_simulated`` column is a simulator-throughput proxy,
    not an operational hardware QUOPS rate.

    The diamond is the rate at the same largest passing shape selected by
    :func:`quops_score_from_summary`.  It is not necessarily the numerically
    largest rate in the table.  Invalid or nonpositive rates are shown as an
    explicit omission count instead of being silently coerced or clipped.
    Main-paper Fig. 2d supplies the decade axes and slope-one reference lines,
    shown here with the notebook's NVIDIA green palette.
    The open circles here retain all measured shape rates, not just score-rate
    pairs from separate architectures. Guides indicate effective executions/s
    after the polarization/acceptance adjustment. Pass ax for a shared figure.
    """
    import matplotlib.pyplot as plt
    required = {"width", "size", rate_column, pass_column}
    if not required.issubset(summary_df.columns):
        raise ValueError(f"summary_df is missing required columns: {sorted(required - set(summary_df.columns))}")
    source_df = summary_df.copy().sort_values(["width", "size"])
    for column in ("width", "size", rate_column):
        source_df[column] = pd.to_numeric(source_df[column], errors="coerce")
    eligible = ~source_df["decision"].eq("invalid") if "decision" in source_df else pd.Series(True, index=source_df.index)
    if "mean_polarization" in source_df:
        eligible &= np.isfinite(pd.to_numeric(source_df["mean_polarization"], errors="coerce"))
    source_df[pass_column] = source_df[pass_column].fillna(False).eq(True) & eligible
    score, best_shape = quops_score_from_summary(source_df, pass_column=pass_column)
    score_symbol = r"\hat{Q}" if "score_passes" in source_df else r"Q_{\mathrm{tested}}"
    valid_mask = np.isfinite(source_df[["width", "size", rate_column]].to_numpy(dtype=float)).all(axis=1) & (source_df[rate_column] > 0.0)
    omitted = int((~valid_mask).sum())
    plot_df = source_df.loc[valid_mask].copy()
    if plot_df.empty:
        raise ValueError(f"{rate_column} does not contain a positive finite rate")
    for row in plot_df[["width", "size"]].itertuples(index=False):
        _positive_integer("width", row.width)
        _positive_integer("size", row.size)
    x_min, x_max = 1.0, 10.0 ** max(1, math.ceil(math.log10(float(source_df["size"].max()))))
    x_padding = 10.0**0.02
    y_min = 10.0 ** math.floor(math.log10(float(plot_df[rate_column].min())))
    y_max = 10.0 ** math.ceil(math.log10(float(plot_df[rate_column].max())))
    y_max = y_min * 10.0 if y_max <= y_min else y_max
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
    _paper_plot_axes(ax)
    ax.set_xlim(x_min / x_padding, x_max * x_padding)
    ax.set_ylim(y_min / x_padding, y_max * x_padding)
    guide_min = math.floor(math.log10(y_min / x_max))
    guide_max = math.ceil(math.log10(y_max / x_min))
    guide_step = max(1, math.ceil((guide_max - guide_min) / 6))
    for exponent in range(guide_min, guide_max + 1, guide_step):
        executions_per_second = 10.0**exponent
        left, right = max(x_min, y_min / executions_per_second), min(x_max, y_max / executions_per_second)
        if right / left <= 1.5:
            continue
        guide_x = np.geomspace(left, right, 80)
        ax.plot(guide_x, executions_per_second * guide_x, color="#B8CCA0", linestyle=(0, (2, 3)), linewidth=0.8, zorder=0)[0].set_gid("quops-rate-guides")
        guide_label_x = left * (right / left)**0.35
        ax.annotate(f"$10^{{{exponent}}}$ eff. circuits/s", (guide_label_x, executions_per_second * guide_label_x), xytext=(0, 4), textcoords="offset points", fontsize=8, color="#888888", ha="center", va="bottom")
    ax.scatter(plot_df["size"], plot_df[rate_column], s=34, marker="o", facecolors="white", edgecolors="#76B900", linewidths=1.1, label="Per-shape rate estimates", zorder=3).set_gid("quops-shape-rates")
    if best_shape is not None:
        selected = plot_df[(plot_df["width"] == best_shape[0]) & (plot_df["size"] == best_shape[1])]
        if not selected.empty:
            selected_rate = float(selected.iloc[0][rate_column])
            ax.scatter([best_shape[1]], [selected_rate], s=80, marker="d", facecolors="#4B7600", edgecolors="#222222", linewidths=0.8, zorder=4, label=f"Rate at ${score_symbol}={score}$: {selected_rate:.3g} QUOPS/sec").set_gid("quops-rate-at-score")
        else:
            ax.text(0.97, 0.04, f"Rate at ${score_symbol}={score}$ unavailable", transform=ax.transAxes, ha="right", fontsize=9, color="#555555")
    else:
        ax.text(0.97, 0.04, "No gated score pass" if "score_passes" in source_df else "No approximate pass in this scan", transform=ax.transAxes, ha="right", fontsize=9, color="#555555")
    ax.set_xlabel("Circuit size, $s$ (QUOPS)", fontsize=11)
    ax.set_ylabel("Rate (QUOPS/sec)", fontsize=11)
    ax.set_title("QUOPS rate by shape" if label is None else f"QUOPS rate by shape\n{label}", fontsize=11)
    if rate_column == "omega_simulated":
        ax.text(0.97, 0.96, "Simulator proxy", transform=ax.transAxes, ha="right", va="top", fontsize=9, color="#555555")
    if omitted:
        ax.text(0.01, 0.01, f"{omitted} invalid/nonpositive rate row(s) omitted", transform=ax.transAxes, fontsize=8, color="#555555", va="bottom")
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.9, fontsize=9)
    if standalone:
        fig.tight_layout()
    return fig


def plot_quops_results(summary_df: pd.DataFrame, *, label: str | None = None):
    """Return the notebook's capability and simulator-rate figure side by side."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), dpi=140)
    if label is not None:
        fig.suptitle(label, fontsize=11)
    _, polarization = plot_quops_scan(summary_df, ax=axes[0])
    plt.close(polarization)
    plot_quops_rate(summary_df, ax=axes[1])
    fig.tight_layout()
    return fig


# ===========================================================================
# Acquisition and provenance: seeds, raw counts and resumable archives
# ===========================================================================

def deterministic_seed(base_seed: int, *coordinates: int) -> int:
    """Derive a coordinate-stable seed without dependence on acquisition order."""
    entropy = [_nonnegative_integer('seed coordinate', value) for value in (base_seed, *coordinates)]
    seed = int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint32)[0])
    return seed if seed else 1


def file_sha256(path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def canonical_json(value) -> str:
    return json.dumps(json_safe(value), sort_keys=True, separators=(',', ':'), allow_nan=False)


def atomic_write_json(path, value) -> None:
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False) + '\n')
    os.replace(temporary, path)


def software_manifest() -> dict:
    versions = {}
    for name in ('numpy', 'pandas', 'matplotlib', 'cudaq', 'cuda-quantum-cu12', 'cuda-quantum-cu13'):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    source_path = Path(__file__)
    return {'python': platform.python_version(), 'platform': platform.platform(), 'cudaq': cudaq.__version__, 'target': str(cudaq.get_target()), 'packages': versions, 'source_sha256': {source_path.name: file_sha256(source_path)}}


def describe_noise(noise, width: int) -> dict:
    """Record static channels on every emitted gate/qubit location in the scan.

    Parameter-dependent callback models need a separately supplied specification;
    this inspection documents the static model used by the bundled examples.
    """
    width = _positive_integer('noise manifest width', width)
    if noise is None:
        return {'model': 'none', 'is_hardware_calibration': False}
    channels = []
    locations = [(gate, [q]) for gate in (*ONE_QUBIT_GATE_KEYS, 'mz') for q in range(width)] + [('cx', [q, r]) for q in range(width) for r in range(width) if q != r]
    for gate, qubits in locations:
        for channel in (noise.get_channels('x', [qubits[1]], [qubits[0]]) if gate == 'cx' else noise.get_channels(gate, qubits)):
            operators = [np.asarray(op, dtype=complex) for op in channel.get_ops()]
            channels.append({'gate': gate, 'qubits': qubits, 'noise_type': str(channel.noise_type), 'parameters': list(channel.parameters), 'kraus': [{'real': op.real.tolist(), 'imag': op.imag.tolist()} for op in operators]})
    return {'model': 'static_cudaq_channels', 'is_hardware_calibration': False, 'channels': channels}


def kernel_manifest(kernel) -> dict:
    gates, q0s, q1s, angles, phis, lams, width = require_kernel_arrays(kernel)
    return {'width': width, 'gates': list(gates), 'q0s': list(q0s), 'q1s': list(q1s), 'angles': list(angles), 'phis': list(phis), 'lams': list(lams)}


def counts_dict(counts) -> dict:
    items = counts.items() if hasattr(counts, 'items') else counts.to_dict().items()
    return {str(bits): _nonnegative_integer('shot count', count) for bits, count in items}


def validate_runtime() -> dict:
    """Check the optimization-boundary feature before collecting any data."""
    try:
        kernel = build_kernel([], [], [], [], 1)
    except TypeError as exc:
        raise RuntimeError('This CUDA-Q build lacks atomic quantum regions. Use the pinned runtime documented in quops.ipynb; do not silently disable MCFE boundaries.') from exc
    if not getattr(kernel, 'atomic_quantum_region', False):
        raise RuntimeError('The CUDA-Q runtime did not preserve the requested atomic quantum region.')
    return {'atomic_quantum_region': True, 'target': str(cudaq.get_target()), 'timing_scope': 'sum_of_M1_M2_sample_call_spans', 'is_hardware_rate': False}


class RunArchive:
    """Archive mirror triplets with paper depth and reject incompatible reuse.

    Shapes and observation rows use the positive even elementary depth ``d``.
    Schema version 2 marks these units explicitly; legacy layer-pair archives
    must be converted separately and are never resumed as paper-depth runs.
    """
    def __init__(self, path, *, shapes, config, noise=None):
        self.path = Path(path)
        self.shapes = default_shape_schedule(shapes)
        for name in ('shots', 'num_circuits', 'num_mirrors'):
            if name in config:
                _positive_integer(name, config[name])
        if not self.shapes:
            raise ValueError('a run requires at least one declared shape')
        if config.get('depth_unit', 'paper') != 'paper':
            raise ValueError("new archives require depth_unit='paper'; convert legacy layer-pair settings explicitly")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.configuration = {'schema_version': 2, 'depth_unit': 'paper', 'config': json_safe(config), 'shapes': self.shapes, 'software': software_manifest(), 'noise': describe_noise(noise, max(shape['width'] for shape in self.shapes))}
        self.fingerprint = sha256(canonical_json(self.configuration).encode()).hexdigest()
        self.config_path = self.path.with_suffix('.config.json')
        self.observations_path = self.path.with_suffix('.observations.jsonl')
        self._records = {}
        if self.config_path.exists():
            saved = json.loads(self.config_path.read_text())
            if saved.get('sha256') != self.fingerprint or canonical_json(saved.get('configuration')) != canonical_json(self.configuration):
                raise ValueError(f'{self.config_path} belongs to different settings, sources, or runtime; select a new run_path')
        elif self.observations_path.exists():
            raise ValueError('orphaned observations have no matching configuration')
        else:
            atomic_write_json(self.config_path, {'sha256': self.fingerprint, 'configuration': self.configuration})
        if self.observations_path.exists():
            for line_number, line in enumerate(self.observations_path.read_text().splitlines(), 1):
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(f'malformed checkpoint line {line_number}; preserve and repair the incomplete line') from exc
                if record.get('config_sha256') != self.fingerprint:
                    raise ValueError('checkpoint configuration does not match the declared run')
                key = self._key(record['row'])
                if key in self._records:
                    raise ValueError(f'duplicate checkpoint {key}')
                self._validate_record(record)
                self._records[key] = record

    @staticmethod
    def _key(row):
        return tuple(_nonnegative_integer(name, row[name]) for name in ('width', 'size', 'circuit_index', 'mirror_index'))

    def _validate_record(self, record):
        row = record['row']
        _paper_shape_to_layer_pairs(row)
        key = self._key(row)
        if (key[0], key[1]) not in {(shape['width'], shape['size']) for shape in self.shapes}:
            raise ValueError('checkpoint shape was not declared')
        shape = next(shape for shape in self.shapes if (shape['width'], shape['size']) == key[:2])
        if _positive_integer('depth', row['depth']) != shape['depth'] or not math.isclose(float(row['zeta']), shape['zeta'], abs_tol=1e-12):
            raise ValueError('checkpoint depth or zeta disagrees with the declared shape')
        config = self.configuration['config']
        if key[2] >= config.get('num_circuits', key[2] + 1) or key[3] >= config.get('num_mirrors', key[3] + 1):
            raise ValueError('checkpoint circuit or mirror index was not declared')
        for ensemble, column in (('br', 'gamma_br'), ('rr', 'gamma_rr'), ('ref', 'gamma_ref')):
            measurement = record['measurements'][ensemble]
            value = effective_polarization_from_counts(measurement['counts'], row['width'], measurement['target'])
            actual = sum(measurement['counts'].values())
            if measurement['shots'] != actual or actual != config.get('shots', actual):
                raise ValueError(f'{ensemble} returned shots differ from the recorded or requested total')
            if not math.isclose(value, float(row[column]), rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError(f'{ensemble} raw counts do not reproduce their saved polarization')
        total = sum(record['measurements'][ensemble]['shots'] for ensemble in ('br', 'rr'))
        if total != _positive_integer('m1_m2_shots', row['m1_m2_shots']):
            raise ValueError('saved M1/M2 shot total disagrees with raw counts')
        if isinstance(row['br_rr_seconds'], bool) or not math.isfinite(float(row['br_rr_seconds'])) or float(row['br_rr_seconds']) <= 0.0:
            raise ValueError('saved M1/M2 sample-call duration must be finite and positive')

    def contains(self, width, size, circuit_index, mirror_index) -> bool:
        return (width, size, circuit_index, mirror_index) in self._records

    def rows(self) -> list[dict]:
        return [record['row'].copy() for record in self._records.values()]

    def record(self, row, *, br_counts, rr_counts, ref_counts, targets, seeds, circuit) -> None:
        measurements = {}
        for ensemble, counts, target in zip(('br', 'rr', 'ref'), (br_counts, rr_counts, ref_counts), targets):
            histogram = counts_dict(counts)
            measurements[ensemble] = {'counts': histogram, 'shots': sum(histogram.values()), 'target': str(target)}
        record = {'config_sha256': self.fingerprint, 'recorded_at_utc': datetime.now(timezone.utc).isoformat(), 'row': json_safe(row), 'seeds': json_safe(seeds), 'source_circuit': kernel_manifest(circuit), 'measurements': measurements}
        self._validate_record(record)
        key = self._key(row)
        if key in self._records:
            raise ValueError(f'checkpoint {key} already exists')
        with self.observations_path.open('a+') as stream:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
            try:
                stream.seek(0)
                for line in stream:
                    if self._key(json.loads(line)['row']) == key:
                        raise ValueError(f'checkpoint {key} was already saved by another writer')
                stream.write(canonical_json(record) + '\n')
                stream.flush()
                os.fsync(stream.fileno())
            finally:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        self._records[key] = record

    def save_summary(self, frame: pd.DataFrame) -> None:
        path = self.path.with_suffix('.summary.csv')
        temporary = path.with_suffix('.csv.tmp')
        frame.to_csv(temporary, index=False)
        os.replace(temporary, path)


# ===========================================================================
# Physical scan runner
# ===========================================================================

def run_quops(*, output='results/quops_run.csv', shapes=None, target='qpp-cpu', shots=1000, num_circuits=50, num_mirrors=1, bootstrap_resamples=2000, seed=2026, p_1q=0.0, p_2q=0.0, familywise_alpha=0.05, estimator='ratio_of_means'):
    """Collect all three MCFE ensembles; resume only an identical declared run.

    Sample times are simulator-call spans, not an operational hardware rate.
    Every shape uses the paper's positive even depth ``d``: ``d // 2`` pairs
    of CNOT and rotation layers. Raw rows, summaries and checkpoints retain
    those units. Shapes from ``default_shape_schedule`` can be passed directly.
    The random seeds and notebook-style acquisition hierarchy are explicit in
    the saved configuration, raw histograms, and original circuit manifests.
    """
    shots = _positive_integer('shots', shots)
    num_circuits = _positive_integer('num_circuits', num_circuits)
    num_mirrors = _positive_integer('num_mirrors', num_mirrors)
    bootstrap_resamples = _positive_integer('bootstrap_resamples', bootstrap_resamples)
    estimator = _validate_estimator(estimator)
    seed = _nonnegative_integer('seed', seed)
    p_1q, p_2q = _noise_probability('p_1q', p_1q), _noise_probability('p_2q', p_2q)
    if (p_1q or p_2q) and target == 'qpp-cpu':
        raise ValueError('Use density-matrix-cpu or a supported NVIDIA target for the noisy example.')
    shapes = default_shape_schedule(shapes if shapes is not None else [{'width': w, 'depth': d, 'zeta': 1.0} for w in range(3, 6) for d in range(4, 12, 2) if in_utility_cone(w, shape_size(w, d, 1.0))])
    validate_inference_config(len(shapes), familywise_alpha=familywise_alpha, bootstrap_resamples=bootstrap_resamples, num_circuits=num_circuits, num_mirrors=num_mirrors)
    previous_target = cudaq.get_target()
    try:
        cudaq.set_target(target)
        validate_runtime()
        noise = make_depolarizing_noise(p_1q, p_2q) if p_1q or p_2q else None
        archive = RunArchive(output, shapes=shapes, config={'target': target, 'shots': shots, 'num_circuits': num_circuits, 'num_mirrors': num_mirrors, 'bootstrap_resamples': bootstrap_resamples, 'seed': seed, 'rc_style': 'pauli_frame', 'confidence_method': 'cluster_bootstrap_normal', 'familywise_alpha': familywise_alpha, 'estimator': estimator}, noise=noise)
        rows = archive.rows()
        for shape in shapes:
            width, depth, size, zeta = shape['width'], shape['depth'], shape['size'], shape['zeta']
            for circuit_index in range(num_circuits):
                circuit_seed = deterministic_seed(seed, width, size, circuit_index)
                circuit = sample_quops_circuit(width, depth, zeta, seed=circuit_seed)
                inverse = adjoint_C(circuit)
                for mirror_index in range(num_mirrors):
                    if archive.contains(width, size, circuit_index, mirror_index):
                        continue
                    rng = np.random.default_rng(deterministic_seed(seed, width, size, circuit_index, mirror_index))
                    cap_seeds, rc_seeds = rng.integers(0, 2**31 - 1, size=3), rng.integers(0, 2**31 - 1, size=3)
                    sample_seed = int(rng.integers(1, 2**31 - 1))
                    caps = [sample_mcfe_caps(width, seed=int(value)) for value in cap_seeds]
                    br_inverse = randomized_compilation(inverse, seed=int(rc_seeds[0]))
                    rr_forward = randomized_compilation(circuit, seed=int(rc_seeds[1]))
                    rr_inverse = randomized_compilation(inverse, seed=int(rc_seeds[2]))
                    cudaq.set_random_seed(sample_seed)
                    started = perf_counter()
                    br_counts = cudaq.sample(BR_kernel, width, caps[0].initial, circuit, br_inverse, caps[0].final, shots_count=shots, noise_model=noise)
                    rr_counts = cudaq.sample(RR_kernel, width, caps[1].initial, rr_forward, rr_inverse, caps[1].final, shots_count=shots, noise_model=noise)
                    elapsed = perf_counter() - started
                    ref_counts = cudaq.sample(REF_kernel, width, caps[2].initial, caps[2].final, shots_count=shots, noise_model=noise)
                    values = [effective_polarization_from_counts(counts, width, cap.target) for counts, cap in zip((br_counts, rr_counts, ref_counts), caps)]
                    row = {'shape': f'w={width},s={size}', 'width': width, 'depth': depth, 'depth_unit': 'paper', 'paper_depth': depth, 'benchmark_depth': depth // 2, 'size': size, 'zeta': zeta, 'circuit_index': circuit_index, 'mirror_index': mirror_index, 'gamma_br': values[0], 'gamma_rr': values[1], 'gamma_ref': values[2], 'gamma': normalized_mcfe_polarization(*values), 'br_rr_seconds': elapsed, 'm1_m2_shots': br_counts.get_total_shots() + rr_counts.get_total_shots()}
                    archive.record(row, br_counts=br_counts, rr_counts=rr_counts, ref_counts=ref_counts, targets=tuple(cap.target for cap in caps), seeds={'circuit': circuit_seed, 'caps': cap_seeds.tolist(), 'rc': rc_seeds.tolist(), 'sample': sample_seed}, circuit=circuit)
                    rows.append(row)
            print(f'completed w={width},s={size}', flush=True)
        data = pd.DataFrame(rows, columns=['shape', 'width', 'depth', 'depth_unit', 'paper_depth', 'benchmark_depth', 'size', 'zeta', 'circuit_index', 'mirror_index', 'gamma_br', 'gamma_rr', 'gamma_ref', 'gamma', 'br_rr_seconds', 'm1_m2_shots'])
        data.to_csv(output, index=False)
        summary = summarize_quops_runs(data, familywise_alpha=familywise_alpha, bootstrap_resamples=bootstrap_resamples, seed=seed, estimator=estimator, declared_shapes=shapes)
        archive.save_summary(summary)
        return data, summary
    finally:
        cudaq.set_target(previous_target)

# ===========================================================================
# Logical plotting and metric displays
# ===========================================================================

def _show(metrics, table, *, title, x, y, xlabel, ax=None, color=None):
    import matplotlib.pyplot as plt
    from IPython.display import display

    display(pd.Series(metrics, name=title))
    if not table.empty:
        table.plot.barh(x=x, y=y, logx=True, legend=False, title=title, xlabel=xlabel, ylabel="", ax=ax, figsize=(7, 3.2) if ax is None else None, color=color)
        if ax is None:
            plt.tight_layout()
            plt.show()


def show_p0(profile, *, ax=None, title="P0: synthesized BR", color="#76B900"):
    """Show a logical profile and its serial depth bound, optionally in a subplot."""

    rows = [{"kind": kind, "operation": operation.split("_", 2)[-1], "count": count} for kind, counts in (("action", profile.actions), ("instrument", profile.instruments)) for operation, count in sorted(counts.items()) if count]
    if profile.idle_sites:
        rows.append({"kind": "memory", "operation": "idle", "count": profile.idle_sites})
    table = pd.DataFrame(rows, columns=["kind", "operation", "count"])
    t_count = sum(count for operation, count in profile.actions.items() if operation.split("_", 2)[-1] in ("t", "tdg"))
    metrics = {"Peak logical qubits": profile.logical_qubits_peak, "T and T-dagger operations": t_count, "Total logical operations": profile.total_operations, "Serial action-depth upper bound": profile.action_depth_upper_bound}
    return _show(metrics, table, title=title, x="operation", y="count", xlabel="Logical operation occurrences (log scale)", ax=ax, color=color)


def plot_p0_comparison(source_resources, resources):
    """Display both P0 summaries and return their before/after synthesis plots."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(20, 4.8), dpi=140, sharex=True)
    show_p0(source_resources, ax=axes[0], title="P0: before synthesis", color="#4B7600")
    show_p0(resources, ax=axes[1], title="P0: after Clifford+T synthesis", color="#76B900")
    fig.tight_layout()
    return fig


def plot_p1_capacity(capacity: pd.DataFrame):
    """Return the compute-region plot using the supplied occupied/free slots."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 3.5), dpi=140)
    capacity.plot.barh(x="region", y=["occupied", "free"], stacked=True, color=["#76B900", "#C9E5A1"], ax=ax, title="Compute region", xlabel="Logical slots", ylabel="")
    fig.tight_layout()
    return fig


def plot_qec_overhead(qec_metrics: pd.DataFrame):
    """Return bars for the checkpoint count and added static work of the BR."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), dpi=140)
    qec_metrics.plot.bar(y="syndrome_extractions", color="#76B900", legend=False, title="Syndrome extraction inside BR", ylabel="Block extractions", rot=0, ax=axes[0])
    qec_metrics.plot.bar(y="additional_instructions", color="#4B7600", legend=False, title="Additional BR instruction cost", ylabel="Static instruction occurrences", rot=0, ax=axes[1])
    for ax in axes:
        ax.set_xlabel("Rounds after each T or T-dagger")
    fig.tight_layout()
    return fig


def show_p2(counts, code):
    """Show static instruction counts and uniform-code block-template costs.

    Every live block is assumed to use the supplied code. External magic-state
    production is excluded from the carrier estimate. Static instruction
    counts include conditional branches and are not primitive gate totals.
    """

    structural = {"call", "repeat", "map_children", "relocate", "establish_support", "establish_topological_record"}
    table = pd.DataFrame([{"instruction": operation, "count": count} for operation, count in sorted(counts.operation_counts.items()) if operation not in structural and count], columns=["instruction", "count"])
    metrics = {"T-state requests": counts.resource_requests.get("t_state", 0), "Syndrome extractions": counts.syndrome_rounds, "Peak encoded blocks": counts.patches_peak, "Carriers in peak block templates": counts.patches_peak * code.block.size}
    return _show(metrics, table, title=f"P2: {code.name} encoded resources", x="instruction", y="count", xlabel="Static Fabric instruction instances (log scale)", color="#76B900")


# ===========================================================================
# Inference calibration with known scalar channels
# ===========================================================================

def _sample_polarizations(rng, polarizations, *, width=2, shots=256):
    values = []
    uniform_shells = np.asarray([math.comb(width, distance) / 2**width for distance in range(width + 1)])
    weights = np.asarray([(-0.5) ** distance for distance in range(width + 1)])
    for polarization in np.asarray(polarizations).flat:
        probabilities = (1.0 - polarization) * uniform_shells
        probabilities[0] += polarization
        counts = rng.multinomial(shots, probabilities)
        moment = float(np.dot(counts, weights) / shots)
        values.append((moment - 4.0**(-width)) / (1.0 - 4.0**(-width)))
    return np.asarray(values).reshape(np.shape(polarizations))


def sample_known_channel_data(rng, truth, *, circuits=60, mirrors=1, shots=256):
    """Use gamma(C)=truth +/-0.12 with exact factorized depolarizing means."""
    target_gamma = truth + rng.choice((-0.12, 0.12), size=circuits)
    reference_gamma = 0.85
    cap_gamma = 0.97
    m1 = np.broadcast_to(target_gamma[:, None] * math.sqrt(reference_gamma * cap_gamma), (circuits, mirrors))
    m2 = np.full((circuits, mirrors), reference_gamma)
    m3 = np.full(circuits * mirrors, cap_gamma)
    br = _sample_polarizations(rng, m1, shots=shots)
    rr = _sample_polarizations(rng, m2, shots=shots)
    ref = _sample_polarizations(rng, m3, shots=shots)
    return ShapePolarizationData(2, 4, 2, 1.0, tuple(map(tuple, br)), tuple(map(tuple, rr)), tuple(ref))


def _wilson(successes, trials):
    z = 1.959963984540054
    rate = successes / trials
    denominator = 1.0 + z * z / trials
    center = (rate + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt(rate * (1.0 - rate) / trials + z * z / (4.0 * trials**2)) / denominator
    return [max(0.0, center - radius), min(1.0, center + radius)]


def calibrate(*, experiments=200, resamples=1000, circuits=60, mirrors=1, shots=256, family_size=5, seed=20260909):
    """Report repeated shape tests, family tests and the prior degenerate case."""
    experiments = _positive_integer("experiments", experiments)
    shots = _positive_integer("shots", shots)
    validate_inference_config(family_size, bootstrap_resamples=resamples, interval_method="percentile", num_circuits=circuits, num_mirrors=mirrors, num_m3_mirrors=circuits * mirrors)
    rng = np.random.default_rng(seed)
    scenarios = (("near_null", QUOPS_THRESHOLD - 0.001, 1), ("below_threshold", 0.45, 1), ("above_threshold", 0.75, 1), ("family_near_null", QUOPS_THRESHOLD - 0.001, family_size))
    results = []
    for scenario, truth, count in scenarios:
        passes = {method: 0 for method in ("normal", "percentile")}
        for _experiment in range(experiments):
            family_passes = {method: False for method in passes}
            for _shape in range(count):
                data = sample_known_channel_data(rng, truth, circuits=circuits, mirrors=mirrors, shots=shots)
                bootstrap_seed = int(rng.integers(0, 2**31 - 1))
                for method in passes:
                    stats = pass_statistics(data, alpha=0.05 / count, resamples=resamples, seed=bootstrap_seed, interval_method=method)
                    family_passes[method] = family_passes[method] or stats["passes"]
            for method in passes:
                passes[method] += int(family_passes[method])
        for method, value in passes.items():
            results.append({"scenario": scenario, "interval_method": method, "true_mean_process_polarization": truth, "experiments": experiments, "family_size": count, "num_circuits": circuits, "num_mirrors": mirrors, "shots_per_mirror": shots, "passes": value, "pass_rate": value / experiments, "monte_carlo_wilson_95": _wilson(value, experiments), "interpretation": "power" if truth > QUOPS_THRESHOLD else "false_pass_probability"})
    # Enumerating the two Bernoulli circuit observations gives an exact check:
    # P(all good)=0.3**2=0.09 formerly gave a false pass with true gamma=0.3.
    exact_false_pass = {method: 0.0 for method in ("normal", "percentile")}
    for first in (0.0, 1.0):
        for second in (0.0, 1.0):
            data = ShapePolarizationData(2, 4, 2, 1.0, ((first, first), (second, second)), ((1.0, 1.0), (1.0, 1.0)), (1.0, 1.0))
            probability = (0.3 if first else 0.7) * (0.3 if second else 0.7)
            for method in exact_false_pass:
                exact_false_pass[method] += probability * pass_statistics(data, resamples=resamples, seed=seed, interval_method=method)["passes"]
    sources = {"src.py": file_sha256(Path(__file__).with_name("src.py"))}
    return {"schema_version": 1, "seed": seed, "bootstrap_resamples": resamples, "alpha_family": 0.05, "threshold": QUOPS_THRESHOLD, "model": "exact depolarizing-channel Hamming-shell sampling; independent target-circuit gamma=truth +/-0.12; M2=0.85; M3=0.97; M1=gamma*sqrt(M2*M3)", "scope": "calibration evidence only for stated scalar model; no hardware, synthesis, MCFE-bias or universal coverage claim", "source_sha256": sources, "results": results, "two_circuit_rare_good_case": {"true_mean_process_polarization": 0.3, "old_degenerate_false_pass_probability": 0.09, "current_exact_false_pass_probability": exact_false_pass}}

# ===========================================================================
# Command-line entry points
# ===========================================================================

def run_quops_cli():
    parser = argparse.ArgumentParser(description="Run physical QUOPS acquisition with durable checkpoints.")
    parser.add_argument('--output', default='results/quops_run.csv')
    parser.add_argument('--target', default='qpp-cpu')
    parser.add_argument('--widths', default='3,4,5')
    parser.add_argument('--depths', default='4,6,8,10', help="Paper depth d: a positive even number of elementary layers, containing d/2 CNOT/rotation pairs")
    parser.add_argument('--zeta', type=float, default=1.0)
    parser.add_argument('--shots', type=int, default=1000)
    parser.add_argument('--circuits', type=int, default=50)
    parser.add_argument('--mirrors', type=int, default=1)
    parser.add_argument('--bootstrap-resamples', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--p-1q', type=float, default=0.0)
    parser.add_argument('--p-2q', type=float, default=0.0)
    parser.add_argument('--familywise-alpha', type=float, default=0.05)
    parser.add_argument('--estimator', choices=MCFE_ESTIMATORS, default='ratio_of_means')
    args = parser.parse_args()
    shapes = [{'width': w, 'depth': d, 'zeta': args.zeta} for w in map(int, args.widths.split(',')) for d in map(int, args.depths.split(',')) if in_utility_cone(w, shape_size(w, d, args.zeta))]
    _, summary = run_quops(output=args.output, shapes=shapes, target=args.target, shots=args.shots, num_circuits=args.circuits, num_mirrors=args.mirrors, bootstrap_resamples=args.bootstrap_resamples, seed=args.seed, p_1q=args.p_1q, p_2q=args.p_2q, familywise_alpha=args.familywise_alpha, estimator=args.estimator)
    import matplotlib.pyplot as plt
    scan, polarization = plot_quops_scan(summary, label=args.target)
    rate = plot_quops_rate(summary, label=args.target)
    for figure, name in ((scan, 'scan'), (polarization, 'polarization'), (rate, 'rate')):
        figure.savefig(Path(args.output).with_suffix(f'.{name}.png'), bbox_inches='tight')
        plt.close(figure)
    print(summary[['width', 'size', 'mean_polarization', 'bootstrap_lower', 'decision']].to_string(index=False))

def calibration_main():
    parser = argparse.ArgumentParser(description="Calibrate QUOPS inference using known scalar channels.")
    parser.add_argument("--experiments", type=int, default=200)
    parser.add_argument("--resamples", type=int, default=1000)
    parser.add_argument("--circuits", type=int, default=60)
    parser.add_argument("--mirrors", type=int, default=1)
    parser.add_argument("--shots", type=int, default=256)
    parser.add_argument("--family-size", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--output", type=Path, default=Path("inference_calibration.json"))
    args = parser.parse_args()
    report = calibrate(experiments=args.experiments, resamples=args.resamples, circuits=args.circuits, mirrors=args.mirrors, shots=args.shots, family_size=args.family_size, seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "results": report["results"], "degenerate_case": report["two_circuit_rare_good_case"]}, indent=2))


def main():
    """Dispatch physical acquisition or known-channel inference calibration."""
    if sys.argv[1:2] == ["calibrate"]:
        del sys.argv[1]
        calibration_main()
    else:
        run_quops_cli()


# ===========================================================================
# Public notebook imports
# ===========================================================================

__all__ = ['BR_kernel', 'DEFAULT_FAMILYWISE_ALPHA', 'MCFE_ESTIMATORS', 'MIN_BOOTSTRAP_TAIL_DRAWS', 'QUOPS_THRESHOLD', 'REF_kernel', 'RR_kernel', 'RunArchive', 'ShapePolarizationData', 'adjoint_C', 'atomic_write_json', 'bootstrap_lower_bound', 'bootstrap_polarization', 'calibrate', 'canonical_json', 'counts_dict', 'default_shape_schedule', 'describe_noise', 'deterministic_seed', 'effective_polarization_from_counts', 'enumerate_shapes', 'file_sha256', 'hamming_distance', 'in_utility_cone', 'json_safe', 'kernel_manifest', 'make_depolarizing_noise', 'make_helios1_noise', 'mcfe_polarization', 'normalized_mcfe_polarization', 'omega_simulated', 'pass_statistics', 'plot_p0_comparison', 'plot_p1_capacity', 'plot_qec_overhead', 'plot_quops_rate', 'plot_quops_results', 'plot_quops_scan', 'plot_quops_timing', 'quops_rate', 'quops_score_from_summary', 'randomized_compilation', 'run_quops', 'sample_known_channel_data', 'sample_mcfe_caps', 'sample_quops_arrays', 'sample_quops_circuit', 'shape_size', 'show_p0', 'show_p2', 'software_manifest', 'summarize_quops_runs', 'utility_cone_bounds', 'validate_inference_config', 'validate_runtime']


if __name__ == "__main__":
    main()
