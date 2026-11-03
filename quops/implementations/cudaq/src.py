"""CUDA-Q implementation of QUOPS for the quops.ipynb tutorial.

Sections: constants; helper functions; circuit construction; randomized
compilation; MCFE caps; noise models; statistical inference; scores and rates;
scan scheduling; plotting; acquisition evidence and archives; logical
resource displays.

The notebook keeps its CUDA-Q demo definitions visible and imports the public
supporting functions listed in __all__.
Logical programs stay visible in the notebook; the helpers here display their
already-compiled resource counts without requiring Logical at module import.

QUOPS estimates mean process polarization over a random circuit ensemble.
Its mirror estimator relies on the reference-compiler assumptions discussed
in the paper. The implemented ratio of ensemble means is approximate;
sampling confidence bounds do not remove its model or normalization bias.
The tutorial reports bootstrap diagnostics and finite-sample component-mean
tests with Hochberg adjustment. Normal-approximation tests remain explicit.
The optional gated-score analysis uses a separate Hochberg region family.
Its two default 5% budgets give a combined 90% region confidence bound.

Use quops.ipynb for QUOPS acquisition and analysis.
"""

from __future__ import annotations

import fcntl
import json
import math
import os
import platform
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata
from pathlib import Path
from typing import Callable, Iterable, Sequence
from uuid import uuid4

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
# Finite-sample bounds cover sampling uncertainty in the MCFE ratio of means.
# The explicit normal approximation and MCFE model assumptions remain separate.
# ===========================================================================

QUOPS_THRESHOLD = 1.0 / math.sqrt(math.e)
MCFE_MIN_DENOMINATOR = 1.0e-12


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


def normalized_mcfe_polarization(gamma_br: float, gamma_rr: float, gamma_ref: float, *, min_denominator: float = MCFE_MIN_DENOMINATOR) -> float:
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
    result = float(gamma_br / (math.sqrt(gamma_rr) * math.sqrt(gamma_ref)))
    return result if math.isfinite(result) else float("nan")


def _resample_mcfe_components(circuit_means, ref, rng, *, resamples):
    """Draw paired BR/RR clusters, then independently resample REF means."""
    num_circuits, num_ref = len(circuit_means), len(ref)
    circuit_indices = rng.integers(0, num_circuits, size=(resamples, num_circuits))
    ref_indices = rng.integers(0, num_ref, size=(resamples, num_ref))
    return circuit_means[circuit_indices], ref[ref_indices].mean(axis=-1)


def _validate_threshold(threshold: float) -> float:
    if isinstance(threshold, (bool, np.bool_)):
        raise ValueError("threshold must be numeric, not boolean")
    threshold = float(threshold)
    if not math.isfinite(threshold) or not 0.0 < threshold <= 1.0:
        raise ValueError("threshold must be finite and in (0, 1]")
    return threshold


# ===========================================================================
# Notebook DataFrame inference
# ===========================================================================

def _normal_p_values(estimates, sigmas, threshold):
    """Evaluate valid one-sided normal tests without cancellation in small tails."""
    with np.errstate(over="ignore"):
        z = (np.asarray(estimates, dtype=float) - np.asarray(threshold, dtype=float)) / np.asarray(sigmas, dtype=float)
    return np.asarray([0.5 * math.erfc(value / math.sqrt(2.0)) for value in z])


def _hochberg_adjusted_p_values(p_values: pd.Series) -> pd.Series:
    """Adjust a complete family's p-values, preserving missing hypotheses as p=1."""
    ordered = p_values.fillna(1.0).sort_values(kind="stable")
    scaled = ordered.to_numpy(dtype=float) * np.arange(len(ordered), 0, -1)
    adjusted = np.minimum(1.0, np.minimum.accumulate(scaled[::-1])[::-1])
    return pd.Series(adjusted, index=ordered.index).reindex(p_values.index)


def _notebook_inference_table(gamma_for_each_shape: pd.DataFrame, sigma_for_each_shape: pd.DataFrame, method: str) -> pd.DataFrame:
    """Match the same observed shapes for diagnostic and certified inference."""
    if method not in ("hoeffding", "normal"):
        raise ValueError("method must be 'hoeffding' or 'normal'")
    _validate_test_shapes(gamma_for_each_shape, "gamma_for_each_shape")
    _validate_test_shapes(sigma_for_each_shape, "sigma_for_each_shape")
    sigma_columns = [column for column in ("width", "size", "sigma", "sigma_valid", "sigma_reason", "bootstrap_valid_fraction", "num_circuits", "num_ref") if column in sigma_for_each_shape]
    result = gamma_for_each_shape.merge(sigma_for_each_shape[sigma_columns], on=["width", "size"], how="outer", validate="one_to_one")
    result["threshold"] = QUOPS_THRESHOLD
    result["inference_method"] = method
    if method == "hoeffding":
        required = {"mean_lambda_br", "mean_lambda_rr", "mean_lambda_ref", "num_circuits", "num_ref"}
        if not required.issubset(result.columns):
            raise ValueError(f"Hoeffding inference requires component means and independent sample counts: missing {sorted(required - set(result.columns))}")
        inverse_dimension_sq = np.power(4.0, -result["width"].to_numpy(dtype=float))
        support_lower = (-0.5 - inverse_dimension_sq) / (1.0 - inverse_dimension_sq)
        means = result[["mean_lambda_br", "mean_lambda_rr", "mean_lambda_ref"]].to_numpy(dtype=float)
        counts = result[["num_circuits", "num_ref"]].to_numpy(dtype=float)
        valid = np.isfinite(means).all(axis=1) & (means >= support_lower[:, None]).all(axis=1) & (means <= 1.0).all(axis=1)
        valid &= np.isfinite(counts).all(axis=1) & (counts >= 1).all(axis=1) & (counts == np.floor(counts)).all(axis=1)
        valid &= (means[:, 1] > MCFE_MIN_DENOMINATOR) & (means[:, 2] > MCFE_MIN_DENOMINATOR)
        result["test_valid"] = valid
        result["test_reason"] = np.where(valid, "finite-sample component-mean bound; assumes independent circuit clusters and independent REF observations", "component-mean bound undefined: missing/out-of-support means, missing sample counts or unstable RR/REF means")
    else:
        result["test_valid"] = np.isfinite(result[["gamma_hat_ws", "sigma"]]).all(axis=1) & result["sigma"].gt(0.0)
        if "sigma_valid" in result:
            result["test_valid"] &= result["sigma_valid"].eq(True)
        result["test_reason"] = "normal one-sided approximation"
        result.loc[~result["test_valid"], "test_reason"] = "normal test undefined: missing estimate or nonpositive/nonfinite uncertainty"
        if "sigma_reason" in result:
            diagnostic = ~result["test_valid"] & result["sigma_reason"].notna()
            result.loc[diagnostic, "test_reason"] = result.loc[diagnostic, "sigma_reason"]
        result.loc[~np.isfinite(result["gamma_hat_ws"]), "test_reason"] = "point estimate undefined: missing data or unstable RR/REF normalization"
    return result


def _hoeffding_mcfe_lower(width: int, mean_br: float, mean_rr: float, mean_ref: float, num_circuits: int, num_ref: int, alpha: float) -> float:
    """Bound the ratio of population means, separately from MCFE model bias.

    Each independent circuit contributes one bounded BR/RR cluster mean. REF
    observations are independent samples. The one-sided error budget is split
    equally across the numerator lower bound and the two denominator upper
    bounds; BR/RR dependence within a circuit does not affect this union bound.
    Positive population RR/REF means are an assumption. A nonpositive numerator
    lower bound returns -infinity, since zero need not bound a negative ratio.
    """
    span = 1.5 / (1.0 - 4.0 ** (-width))

    # Add a small outward margin for floating-point evaluation of the radius.
    rounding_margin = 1.0 + 8.0 * np.finfo(float).eps
    radius_circuit = math.nextafter(rounding_margin * span * math.sqrt(math.log(3.0 / alpha) / (2.0 * num_circuits)), math.inf)
    radius_ref = math.nextafter(rounding_margin * span * math.sqrt(math.log(3.0 / alpha) / (2.0 * num_ref)), math.inf)
    lower_br = math.nextafter(mean_br - radius_circuit, -math.inf)
    upper_rr = min(1.0, math.nextafter(mean_rr + radius_circuit, math.inf))
    upper_ref = min(1.0, math.nextafter(mean_ref + radius_ref, math.inf))
    if lower_br <= 0.0 or upper_rr <= 0.0 or upper_ref <= 0.0:
        return -math.inf
    return math.nextafter(lower_br / (math.sqrt(upper_rr) * math.sqrt(upper_ref)), -math.inf)


def _hoeffding_mcfe_p_value(width: int, mean_br: float, mean_rr: float, mean_ref: float, num_circuits: int, num_ref: int, threshold: float) -> float:
    """Invert the monotone bound and retain the upper crossing endpoint."""
    arguments = (width, mean_br, mean_rr, mean_ref, num_circuits, num_ref)
    if _hoeffding_mcfe_lower(*arguments, 1.0) <= threshold:
        return 1.0
    lower, upper = 0.0, 1.0
    for _ in range(70):
        midpoint = (lower + upper) / 2.0
        if _hoeffding_mcfe_lower(*arguments, midpoint) > threshold:
            upper = midpoint
        else:
            lower = midpoint
    return min(1.0, math.nextafter(upper, 1.0))


def calculate_gamma_lower(gamma_for_each_shape: pd.DataFrame, sigma_for_each_shape: pd.DataFrame, *, method: str = "hoeffding", alpha: float = 0.05) -> pd.DataFrame:
    """Return a one-sided bound and classify individual shape tests.

    The default Hoeffding bound has finite-sample coverage for the ratio of
    population BR/RR/REF means under independent circuit clusters, independent
    REF observations and positive population denominator means. It does not
    remove the physical/reference-compiler approximation in MCFE. The bootstrap
    sigma remains a diagnostic. method='normal' selects the paper's normal
    approximation explicitly; its finite-sample coverage is not guaranteed.
    """
    alpha = _probability("alpha", alpha)
    result = _notebook_inference_table(gamma_for_each_shape, sigma_for_each_shape, method)
    result["gamma_lower"] = np.nan
    if method == "hoeffding":
        for index in result.index[result["test_valid"]]:
            row = result.loc[index]
            result.loc[index, "gamma_lower"] = _hoeffding_mcfe_lower(int(row["width"]), row["mean_lambda_br"], row["mean_lambda_rr"], row["mean_lambda_ref"], int(row["num_circuits"]), int(row["num_ref"]), alpha)
    else:
        from statistics import NormalDist
        result.loc[result["test_valid"], "gamma_lower"] = result.loc[result["test_valid"], "gamma_hat_ws"] - NormalDist().inv_cdf(1.0 - alpha) * result.loc[result["test_valid"], "sigma"]
    result["interval_valid"] = result["test_valid"]
    result["confidence"] = 1.0 - alpha
    result["passes"] = result["test_valid"] & result["gamma_lower"].gt(result["threshold"])
    result["decision"] = np.where(result["test_valid"], np.where(result["passes"], "Pass", "Fail"), "Invalid")
    return result

def _validate_test_shapes(data: pd.DataFrame, name: str) -> None:
    """Require one row per positive integer (width, size) without silently dropping shapes."""
    if not {"width", "size"}.issubset(data.columns):
        raise ValueError(f"{name} requires width and size columns")
    for column in ("width", "size"):
        data[column].map(lambda value: _positive_integer(column, value))
    if data.duplicated(["width", "size"]).any():
        raise ValueError(f"{name} contains duplicate shapes")

def calculate_p_value(gamma_for_each_shape: pd.DataFrame, sigma_for_each_shape: pd.DataFrame, *, method: str = "hoeffding") -> pd.DataFrame:
    """Return shape p-values for the selected sampling-inference method.

    The default inverts the finite-sample component-mean bound. Its null is
    that the population MCFE ratio is at most the QUOPS threshold; interpreting
    that ratio as physical process polarization still needs the MCFE model.
    method='normal' returns the approximate SI Eq. (61) test. Invalid shapes
    remain present with NaN p-values and cannot pass the family procedure.
    """
    result = _notebook_inference_table(gamma_for_each_shape, sigma_for_each_shape, method)
    result["p_value"] = np.nan
    if method == "hoeffding":
        for index in result.index[result["test_valid"]]:
            row = result.loc[index]
            result.loc[index, "p_value"] = _hoeffding_mcfe_p_value(int(row["width"]), row["mean_lambda_br"], row["mean_lambda_rr"], row["mean_lambda_ref"], int(row["num_circuits"]), int(row["num_ref"]), QUOPS_THRESHOLD)
    else:
        result.loc[result["test_valid"], "p_value"] = _normal_p_values(result.loc[result["test_valid"], "gamma_hat_ws"], result.loc[result["test_valid"], "sigma"], QUOPS_THRESHOLD)
    return result

def hochberg_test(p_value_for_each_shape: pd.DataFrame, alpha: float = 0.05, *, declared_shapes: Iterable[dict] | pd.DataFrame | None = None) -> pd.DataFrame:
    """Apply SI Eq. (67) to one complete family and return every shape's decision.

    Use all planned shapes, chosen independently of the test data. Supplying
    declared_shapes retains missing shapes as Invalid; otherwise the input
    rows define the family. Invalid tests occupy a place with p=1 for the
    adjustment and cannot pass. Do not filter to individual passes first.
    Hochberg's FWER control assumes valid p-values and independence or the
    positive dependence described in the paper. This helper does not gate
    a score sequence or combine separate score and region error budgets.
    """
    alpha = _probability("alpha", alpha)
    _validate_test_shapes(p_value_for_each_shape, "p_value_for_each_shape")
    if "p_value" not in p_value_for_each_shape:
        raise ValueError("p_value_for_each_shape requires a p_value column")
    result = p_value_for_each_shape.copy()
    if declared_shapes is not None:
        declared = pd.DataFrame(declared_shapes)
        _validate_test_shapes(declared, "declared_shapes")
        supplied = set(result[["width", "size"]].itertuples(index=False, name=None))
        planned = set(declared[["width", "size"]].itertuples(index=False, name=None))
        if not supplied.issubset(planned):
            raise ValueError("p-values contain shapes outside declared_shapes")
        result = declared[["width", "size"]].merge(result, on=["width", "size"], how="left", validate="one_to_one")
    if result.empty:
        raise ValueError("Hochberg requires at least one declared shape")
    result = result.reset_index(drop=True)
    p_values = pd.to_numeric(result["p_value"], errors="raise")
    if (p_values.notna() & ~p_values.between(0.0, 1.0)).any():
        raise ValueError("p-values must be in [0, 1] or NaN for invalid tests")
    valid = p_values.notna()
    if "test_valid" in result:
        valid &= result["test_valid"].eq(True)
    result["p_value"] = p_values
    result["test_valid"] = valid
    if "test_reason" in result:
        result["test_reason"] = result["test_reason"].fillna("no observations for this declared shape")

    # Every declared hypothesis participates, including invalid cases with p=1.
    result["adjusted_p_value"] = _hochberg_adjusted_p_values(p_values.where(valid, 1.0))
    result["passes"] = valid & result["adjusted_p_value"].le(alpha)
    result["decision"] = np.where(valid, np.where(result["passes"], "Pass", "Fail"), "Invalid")
    result["alpha"] = alpha
    result["family_size"] = len(result)
    result["multiplicity"] = "hochberg"
    if "threshold" in result:
        result["threshold"] = result["threshold"].fillna(QUOPS_THRESHOLD)
    else:
        result["threshold"] = QUOPS_THRESHOLD
    return result


def _validate_score_sequence(sequence, declared_shapes):
    """Check SI Eqs. (63)--(64); the caller must choose S before its test data."""
    try:
        sequence = [(_positive_integer("score width", width), _positive_integer("score size", size)) for width, size in sequence]
    except (TypeError, ValueError) as exc:
        raise ValueError("score_sequence must contain (width, size) integer pairs") from exc
    declared = set(declared_shapes[["width", "size"]].itertuples(index=False, name=None))
    if not sequence or any(shape not in declared or not in_utility_cone(*shape) for shape in sequence):
        raise ValueError("score_sequence must contain declared shapes inside the utility cone")
    if any(right[1] <= left[1] for left, right in zip(sequence, sequence[1:])):
        raise ValueError("score_sequence sizes must be strictly increasing in their predeclared order")
    return sequence


def gated_hochberg_test(p_value_for_each_shape: pd.DataFrame, score_sequence: Sequence[tuple[int, int]], *, declared_shapes: Iterable[dict] | pd.DataFrame, score_alpha: float = 0.05, region_alpha: float = 0.05) -> pd.DataFrame:
    """Apply SI Eqs. (63)--(68) to the current notebook's shape p-values.

    Fix the declared family and increasing in-cone sequence S independently
    of the final test data. This function checks the declaration's structure;
    it cannot establish when the caller chose it. Each reached shape in S is
    tested at score_alpha. The first failure or invalid test closes the gate.
    All later S shapes remain Not tested, even if observations were collected.

    R contains every declared shape outside S, including missing or invalid
    hypotheses. Only R receives the separate Hochberg adjustment. score_passes
    selects the score; passes combines successful S and R shapes. Its default
    joint error budget is 0.10, not 0.05. Extending these passes to a downward
    closed capability region additionally assumes that mean polarization is
    non-increasing in both width and size.

    The input can use the conservative Hoeffding or explicit normal method.
    Valid shape tests, the paper's dependence conditions, and the MCFE model
    remain assumptions. Normal p-values retain their approximation limits.
    """
    score_alpha = _probability("score_alpha", score_alpha)
    region_alpha = _probability("region_alpha", region_alpha)
    if score_alpha + region_alpha >= 1.0:
        raise ValueError("score_alpha + region_alpha must be less than 1")
    if declared_shapes is None:
        raise ValueError("gated testing requires the complete declared_shapes family")
    declared = pd.DataFrame(declared_shapes)
    _validate_test_shapes(declared, "declared_shapes")
    sequence = _validate_score_sequence(score_sequence, declared)
    threshold_values = pd.to_numeric(p_value_for_each_shape["threshold"], errors="raise").dropna().unique() if "threshold" in p_value_for_each_shape else np.array([])
    if len(threshold_values) > 1:
        raise ValueError("gated testing requires one common threshold")
    threshold = _validate_threshold(threshold_values[0]) if len(threshold_values) else QUOPS_THRESHOLD
    methods = p_value_for_each_shape["inference_method"].dropna().unique() if "inference_method" in p_value_for_each_shape else np.array([])
    if len(methods) > 1:
        raise ValueError("gated testing requires one sampling-inference method")
    if len(methods) and methods[0] not in ("hoeffding", "normal"):
        raise ValueError("inference_method must be 'hoeffding' or 'normal'")

    # Reuse the existing complete-family validation, then apply S and R separately.
    result = hochberg_test(p_value_for_each_shape, alpha=region_alpha, declared_shapes=declared)
    result["threshold"] = threshold
    if len(methods):
        result["inference_method"] = methods[0]
    elif "inference_method" in result:
        result = result.drop(columns=["inference_method"])
    if "test_reason" not in result:
        result["test_reason"] = np.where(result["test_valid"], "supplied shape p-value", "missing or invalid shape p-value")
    indices = {(int(row["width"]), int(row["size"])): index for index, row in result.iterrows()}
    result["adjusted_p_value"] = np.nan
    result["test_family"] = "Hochberg (R)"
    result["tested"] = True
    result["score_passes"] = False
    result["passes"] = False
    result["decision"] = "Invalid"
    result["decision_reason"] = "missing or invalid shape test"
    result["declared_family_size"] = len(result)
    result["family_size"] = len(result) - len(sequence)

    # The first non-pass stops S; its unused tail never enters R.
    gate_open = True
    for shape in sequence:
        index = indices[shape]
        result.loc[index, "test_family"] = "Gated score (S)"
        result.loc[index, "alpha"] = score_alpha
        result.loc[index, "family_size"] = len(sequence)
        if not gate_open:
            result.loc[index, ["p_value", "adjusted_p_value"]] = np.nan
            result.loc[index, ["test_valid", "tested"]] = False
            result.loc[index, "decision"] = "Not tested"
            result.loc[index, "decision_reason"] = "after the first unsuccessful gated test"
            result.loc[index, "test_reason"] = "not tested after the first unsuccessful gated test"
            continue
        valid = bool(result.loc[index, "test_valid"])
        passed = bool(valid and result.loc[index, "p_value"] <= score_alpha)
        result.loc[index, ["score_passes", "passes"]] = passed
        result.loc[index, "decision"] = "Pass" if passed else "Fail" if valid else "Invalid"
        result.loc[index, "decision_reason"] = "shape p-value compared with the gated score budget" if valid else "missing or invalid shape test stops the score sequence"
        gate_open = passed

    # R keeps all its declared hypotheses, with p=1 for missing or invalid tests.
    remaining = result.index[result["test_family"].eq("Hochberg (R)")]
    if len(remaining):
        valid = result.loc[remaining, "test_valid"]
        result.loc[remaining, "adjusted_p_value"] = _hochberg_adjusted_p_values(result.loc[remaining, "p_value"].where(valid, 1.0))
        result.loc[remaining, "passes"] = valid & result.loc[remaining, "adjusted_p_value"].le(region_alpha)
        result.loc[remaining, "decision"] = np.where(valid, np.where(result.loc[remaining, "passes"], "Pass", "Fail"), "Invalid")
        result.loc[remaining, "decision_reason"] = np.where(valid, "Hochberg-adjusted p-value compared with the region budget", "missing or invalid shape test remains in the region family")
    result["multiplicity"] = "gated_hochberg"
    result["score_alpha"] = score_alpha
    result["region_alpha"] = region_alpha
    result["score_confidence"] = 1.0 - score_alpha
    result["region_confidence"] = 1.0 - score_alpha - region_alpha
    result.attrs.update({"multiplicity": "gated_hochberg", "score_sequence": sequence, "score_alpha": score_alpha, "region_alpha": region_alpha, "score_confidence": 1.0 - score_alpha, "region_confidence": 1.0 - score_alpha - region_alpha, "region_assumption": "downward closure requires mean polarization non-increasing in width and size"})
    return result


def _validate_lambda_observations(data: pd.DataFrame) -> None:
    """Require finite observations and balanced circuit clusters before inference."""
    columns = ["width", "size", "circuit_index", "lambda_br", "lambda_rr", "lambda_ref"]
    missing = set(columns) - set(data.columns)
    if missing:
        raise ValueError(f"data is missing required columns: {sorted(missing)}")
    if data.empty or data[columns].isna().any().any():
        raise ValueError("data must contain nonmissing shape, circuit and lambda observations")
    if not np.isfinite(data[["lambda_br", "lambda_rr", "lambda_ref"]].to_numpy(dtype=float)).all():
        raise ValueError("lambda observations must be finite")
    for column in ("width", "size", "circuit_index"):
        validator = _nonnegative_integer if column == "circuit_index" else _positive_integer
        data[column].map(lambda value: validator(column, value))
    inverse_dimension_sq = np.power(4.0, -data["width"].to_numpy(dtype=float))
    lower = (-0.5 - inverse_dimension_sq) / (1.0 - inverse_dimension_sq)
    values = data[["lambda_br", "lambda_rr", "lambda_ref"]].to_numpy(dtype=float)
    if np.any(values < lower[:, None]) or np.any(values > 1.0):
        raise ValueError("effective polarizations are outside the bounded Hamming-observable support")
    if "mirror_index" in data:
        data["mirror_index"].map(lambda value: _nonnegative_integer("mirror_index", value))
        if data.duplicated(["width", "size", "circuit_index", "mirror_index"]).any():
            raise ValueError("data contains duplicate circuit/mirror observations")
    mirror_counts = data.groupby(["width", "size", "circuit_index"]).size()
    if mirror_counts.groupby(level=["width", "size"]).nunique().gt(1).any():
        raise ValueError("Each shape must have the same number of mirrors for every original circuit; complete the acquisition before analysis.")


def validate_acquisition(data: pd.DataFrame, declared_shapes: Iterable[dict], *, num_circuits: int, num_mirrors: int, shots: int) -> None:
    """Require the complete fixed notebook acquisition before drawing conclusions.

    The declaration must be fixed before observing results. Every planned
    original circuit and mirror index must occur exactly once, with the
    requested number of shots in each of BR, RR and REF. This detects an
    interrupted prefix even when its observed circuit groups are balanced.
    """
    num_circuits = _positive_integer("num_circuits", num_circuits)
    num_mirrors = _positive_integer("num_mirrors", num_mirrors)
    shots = _positive_integer("shots", shots)
    shapes = default_shape_schedule(declared_shapes)
    if not shapes:
        raise ValueError("acquisition requires at least one declared shape")
    required = {"width", "size", "depth", "zeta", "circuit_index", "mirror_index", "shots_br", "shots_rr", "shots_ref", "m1_m2_shots"}
    missing_columns = required - set(data.columns)
    if missing_columns:
        raise ValueError(f"acquisition data is missing columns: {sorted(missing_columns)}")
    _validate_lambda_observations(data)
    expected = {(shape["width"], shape["size"], circuit, mirror) for shape in shapes for circuit in range(num_circuits) for mirror in range(num_mirrors)}
    observed = set(data[["width", "size", "circuit_index", "mirror_index"]].itertuples(index=False, name=None))
    if observed != expected:
        raise ValueError(f"acquisition is incomplete or differs from its declaration: {len(expected - observed)} missing and {len(observed - expected)} unexpected circuit/mirror observations")
    planned = {(shape["width"], shape["size"]): shape for shape in shapes}
    for row in data.to_dict("records"):
        shape = planned[(row["width"], row["size"])]
        _paper_shape_to_layer_pairs(row)
        if _positive_integer("depth", row["depth"]) != shape["depth"] or not math.isclose(float(row["zeta"]), shape["zeta"], rel_tol=0.0, abs_tol=1.0e-12) or shape_size(row["width"], row["depth"], row["zeta"]) != row["size"]:
            raise ValueError("observed depth, final-layer fraction or size differs from the declared shape")
        for ensemble in ("br", "rr", "ref"):
            if _nonnegative_integer(f"shots_{ensemble}", row[f"shots_{ensemble}"]) != shots:
                raise ValueError(f"{ensemble} returned shot count differs from the declared shots")
        if _nonnegative_integer("m1_m2_shots", row["m1_m2_shots"]) != row["shots_br"] + row["shots_rr"]:
            raise ValueError("M1/M2 shot total differs from the individual BR/RR counts")


def calculate_mean_lambdas(data: pd.DataFrame) -> pd.DataFrame:
    """Average finite lambdas for each shape, requiring equal mirror counts per circuit."""
    _validate_lambda_observations(data)
    return data.groupby(["width", "size"], as_index=False).agg(mean_lambda_br=("lambda_br", "mean"), mean_lambda_rr=("lambda_rr", "mean"), mean_lambda_ref=("lambda_ref", "mean"))

def calculate_gamma_hat_ws(mean_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate the estimated mean process polarization for each shape."""
    result = mean_data.copy()
    result["gamma_hat_ws"] = [normalized_mcfe_polarization(br, rr, ref) for br, rr, ref in result[["mean_lambda_br", "mean_lambda_rr", "mean_lambda_ref"]].itertuples(index=False, name=None)]
    return result

def calculate_quops_rate(data: pd.DataFrame, gamma_for_each_shape: pd.DataFrame) -> pd.DataFrame:
    """Return simulator QUOPS/s from each shape's summed BR/RR time and shots.

    Optional m1_m2_kept_shots records retained BR/RR shots per observation.
    Validate those counts before summing, then apply Eq. (72) to the shape's
    total and retained shots. Omitting the column means all shots were kept.
    The times remain simulator sampling-call spans, not hardware acquisition.
    The square-root acceptance factor is the paper's mirror-half heuristic;
    this calculation does not validate the experiment's postselection rule.
    """
    has_kept_shots = "m1_m2_kept_shots" in data
    measurements = data[["width", "size", "br_rr_seconds", "m1_m2_shots"] + (["m1_m2_kept_shots"] if has_kept_shots else [])].copy()
    measurements["m1_m2_shots"] = measurements["m1_m2_shots"].map(lambda value: _nonnegative_integer("m1_m2_shots", value))
    if has_kept_shots:
        measurements["m1_m2_kept_shots"] = measurements["m1_m2_kept_shots"].map(lambda value: _nonnegative_integer("m1_m2_kept_shots", value))
        if measurements["m1_m2_kept_shots"].gt(measurements["m1_m2_shots"]).any():
            raise ValueError("m1_m2_kept_shots cannot exceed total m1_m2_shots in any observation")
    if not np.isfinite(measurements["br_rr_seconds"]).all() or not measurements["br_rr_seconds"].gt(0.0).all():
        raise ValueError("br_rr_seconds must contain finite, positive sampling times")

    # Sum every BR/RR observation before matching its shape's gamma; REF time and shots are excluded.
    aggregations = {"tau_wall_seconds": ("br_rr_seconds", "sum"), "m1_m2_shots": ("m1_m2_shots", "sum")}
    if has_kept_shots:
        aggregations["m1_m2_kept_shots"] = ("m1_m2_kept_shots", "sum")
    totals = measurements.groupby(["width", "size"], as_index=False).agg(**aggregations)
    result = totals.merge(gamma_for_each_shape[["width", "size", "gamma_hat_ws"]], on=["width", "size"], how="left", validate="one_to_one")

    # Use the measured gamma directly, including sampling fluctuations above one; undefined or negative estimates yield NaN.
    result["omega_simulated"] = [omega_simulated(row.size, row.gamma_hat_ws, row.tau_wall_seconds, row.m1_m2_shots, kept_shots=row.m1_m2_kept_shots if has_kept_shots else row.m1_m2_shots) if np.isfinite(row.gamma_hat_ws) and row.gamma_hat_ws >= 0.0 else np.nan for row in result.itertuples(index=False)]
    return result

def calculate_sigma(data: pd.DataFrame, bootstrap_resamples: int = 2000, seed: int = 2026) -> pd.DataFrame:
    """Estimate gamma uncertainty while retaining untestable shapes.

    Point estimates and every bootstrap draw use the same finite/tiny
    denominator rule. If any draw is invalid, sigma is NaN: invalid draws
    are never discarded. Diagnostics distinguish this case from zero
    bootstrap spread and insufficient independent circuit observations. The
    independent circuit and REF counts also support the separate finite-sample
    component-mean test; mirrors never increase the circuit-cluster count.
    """
    _validate_lambda_observations(data)
    bootstrap_resamples = _positive_integer("bootstrap_resamples", bootstrap_resamples)
    if bootstrap_resamples < 2:
        raise ValueError("At least two bootstrap resamples are required to estimate uncertainty.")
    rng = np.random.default_rng(seed)
    rows = []

    for (width, size), group in data.groupby(["width", "size"]):
        circuit_means = group.groupby("circuit_index")[["lambda_br", "lambda_rr"]].mean().to_numpy()
        ref_lambdas = group["lambda_ref"].to_numpy()
        num_circuits = len(circuit_means)
        counts = {"num_circuits": num_circuits, "num_ref": len(ref_lambdas)}
        if num_circuits < 2:
            rows.append({"width": width, "size": size, **counts, "sigma": np.nan, "sigma_valid": False, "sigma_reason": "fewer than two independent original circuits", "bootstrap_valid_fraction": np.nan})
            continue

        # Preserve paired circuit clusters, independent REF draws and the notebook's batch RNG order.
        sampled, bootstrap_ref = _resample_mcfe_components(circuit_means, ref_lambdas, rng, resamples=bootstrap_resamples)
        bootstrap_means = sampled.mean(axis=1)
        bootstrap_br = bootstrap_means[:, 0]
        bootstrap_rr = bootstrap_means[:, 1]

        # Keep every draw in the diagnostic; one unstable normalization makes this shape untestable.
        valid_draws = np.isfinite(bootstrap_br) & np.isfinite(bootstrap_rr) & np.isfinite(bootstrap_ref) & (bootstrap_rr > MCFE_MIN_DENOMINATOR) & (bootstrap_ref > MCFE_MIN_DENOMINATOR)
        bootstrap_gamma = np.full(bootstrap_resamples, np.nan)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            bootstrap_gamma[valid_draws] = bootstrap_br[valid_draws] / (np.sqrt(bootstrap_rr[valid_draws]) * np.sqrt(bootstrap_ref[valid_draws]))
        valid_draws &= np.isfinite(bootstrap_gamma)
        if not valid_draws.all():
            rows.append({"width": width, "size": size, **counts, "sigma": np.nan, "sigma_valid": False, "sigma_reason": "nonfinite or tiny RR/REF normalization in one or more bootstrap draws", "bootstrap_valid_fraction": float(valid_draws.mean())})
            continue
        with np.errstate(over="ignore", invalid="ignore"):
            sigma = float(np.std(bootstrap_gamma, ddof=1))
        valid_sigma = math.isfinite(sigma) and sigma > 0.0
        reason = "valid cluster-bootstrap standard deviation" if valid_sigma else "zero or nonfinite bootstrap spread; the normal test is undefined"
        rows.append({"width": width, "size": size, **counts, "sigma": sigma, "sigma_valid": valid_sigma, "sigma_reason": reason, "bootstrap_valid_fraction": 1.0})

    return pd.DataFrame(rows)


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
    were scheduled and tested. Gated results use score_passes by default;
    successful region-family tests cannot increase the gated score.
    """
    if pass_column == "passes" and "score_passes" in summary_df:
        pass_column = "score_passes"
    required = {"width", "size", pass_column}
    if summary_df.empty or not required.issubset(summary_df.columns):
        return 0, None
    eligible = summary_df[pass_column].fillna(False).eq(True)
    if "decision" in summary_df:
        eligible &= summary_df["decision"].astype(str).str.casefold().eq("pass")
    if "measured" in summary_df:
        eligible &= summary_df["measured"].fillna(False).eq(True)
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

    ``gamma`` is a point estimate. Sampling fluctuations can put it above one;
    Eq. (72) uses that estimate directly, without clipping or discarding it.

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
    if gamma < 0.0:
        raise ValueError("gamma used for a rate must be non-negative")
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
# Scan scheduling
# ===========================================================================

def _default_shape_schedule_from_layer_pairs(shapes: Iterable[dict] | Sequence[int]) -> list[dict]:
    """Validate and order the complete finite family before acquisition.

    This deterministic order defines the complete declared acquisition family.
    Statistical testing uses the separate notebook inference functions.

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
    """Validate public paper-depth metadata and derive internal layer pairs."""
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


def _paper_plot_axes(ax):
    """Apply the main paper's decade scales, fine ticks and uncluttered axes."""
    from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullFormatter
    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10, numticks=8))
        axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
        axis.set_major_formatter(LogFormatterMathtext(base=10))
        axis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", direction="out", length=4, width=0.8, labelsize=10, color="#555555")
    ax.tick_params(which="minor", direction="out", length=2, width=0.6, color="#777777")
    for spine in ax.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(0.8)
    ax.grid(False, which="both")
    ax.set_axisbelow(True)


def _highlight_quops_point(ax, x: float, y: float, text: str, *, gid: str, legend_label: str) -> None:
    """Mark the selected shape and place its numerical label toward the plot interior."""
    ax.scatter([x], [y], s=165, marker="D", facecolors="#4B7600", edgecolors="#203400", linewidths=1.3, label=legend_label, zorder=6).set_gid(gid)
    fraction_x, fraction_y = ax.transAxes.inverted().transform(ax.transData.transform((x, y)))
    offset_x = 26 if fraction_x < 0.45 else -26
    offset_y = 34 if fraction_y < 0.5 else -34
    horizontal = "left" if offset_x > 0 else "right"
    vertical = "bottom" if offset_y > 0 else "top"
    ax.annotate(text, xy=(x, y), xytext=(offset_x, offset_y), textcoords="offset points", ha=horizontal, va=vertical, fontsize=10, fontweight="bold", color="#203400", bbox={"boxstyle": "round,pad=0.45", "facecolor": "#F6FAEF", "edgecolor": "#4B7600", "alpha": 0.98}, arrowprops={"arrowstyle": "->", "color": "#4B7600", "lw": 1.3}, zorder=7).set_gid(f"{gid}-label")


def _plot_quops_individual_scan(summary_df: pd.DataFrame, *, title: str, label: str | None, threshold: float, figsize: tuple[float, float], dpi: int, ax=None, test_description: str = "Individual 95% tests", assume_monotonic: bool = False):
    """Draw supplied shape decisions and the geometric utility cone."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    plot_df = summary_df.copy().sort_values(["width", "size"])
    if plot_df.empty:
        raise ValueError("summary_df does not contain plottable rows")
    for column in ("width", "size"):
        plot_df[column] = plot_df[column].map(lambda value: _positive_integer(column, value))
    plot_df["decision"] = plot_df["decision"].astype(str).str.capitalize()
    if not plot_df["decision"].isin(["Pass", "Fail", "Invalid", "Not tested"]).all():
        raise ValueError("shape decisions must be Pass, Fail, Invalid or Not tested")

    # Shade the geometric utility cone; each marker represents one measured shape.
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, layout="constrained")
    else:
        fig = ax.figure
    size_grid = np.geomspace(1.0 if assume_monotonic else max(1.0, plot_df["size"].min() / 1.8), plot_df["size"].max() * 1.8, 300)
    ax.fill_between(size_grid, size_grid ** (1 / 3), np.sqrt(size_grid), color="#EDF5DD", label=r"Utility cone: $w^2 \leq s \leq w^3$", zorder=0).set_gid("quops-utility-cone")
    ax.plot(size_grid, np.sqrt(size_grid), color="#91AD63", linewidth=1.0, linestyle="--", label=r"$s=w^2$")
    ax.plot(size_grid, size_grid ** (1 / 3), color="#91AD63", linewidth=1.0, linestyle=":", label=r"$s=w^3$")

    # Eq. (68) adds smaller shapes only under the explicit monotonicity assumption.
    if assume_monotonic:
        passes = plot_df.loc[plot_df["decision"].eq("Pass")]
        if not passes.empty:
            region_x = sorted({1.0, float(size_grid[-1]), *passes["size"].astype(float)})
            region_y = [max([1.0, *passes.loc[passes["size"].ge(right), "width"].astype(float)]) for right in region_x[1:]] + [1.0]
            confidence = float(plot_df["region_confidence"].iloc[0]) if plot_df["multiplicity"].iloc[0] == "gated_hochberg" else 1.0 - float(plot_df["alpha"].iloc[0])
            ax.fill_between(region_x, 1.0, region_y, step="post", color="#B4D984", alpha=0.4, label=f"Inferred region ({confidence:.0%}; monotonicity assumed)", zorder=0.5).set_gid("quops-capability-region")

    # Use the supplied test decisions to distinguish the tested shapes.
    for decision, color, marker in [("Pass", "#76B900", "o"), ("Fail", "#C5523F", "X"), ("Invalid", "#777777", "s"), ("Not tested", "#999999", "^")]:
        points = plot_df.loc[plot_df["decision"].eq(decision)]
        if not points.empty:
            marker_label = "Threshold not established" if decision == "Fail" else "Not tested: gate closed" if decision == "Not tested" else decision
            ax.scatter(points["size"], points["width"], s=85, marker=marker, color=color, edgecolors="white", linewidths=0.8, label=marker_label, zorder=3).set_gid(f"quops-{decision.lower().replace(' ', '-')}")

    # Show circuit size horizontally and qubit width vertically, as in the paper.
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(size_grid[0], size_grid[-1])
    ax.set_ylim(min(size_grid[0] ** (1 / 3), plot_df["width"].min()) / 1.15, max(np.sqrt(size_grid[-1]), plot_df["width"].max()) * 1.15)
    ax.set_xticks(sorted(plot_df["size"].unique()))
    ax.set_yticks(sorted(plot_df["width"].unique()))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.minorticks_off()
    ax.set_xlabel("Circuit size, s (QUOPS)")
    ax.set_ylabel("Circuit width, w (qubits)")
    heading = title if label is None else f"{title}\n{label}"
    ax.set_title(f"{heading}\n{test_description}; threshold = {threshold:.4f}")
    ax.grid(alpha=0.2)
    ax.set_axisbelow(True)

    # Select the largest passing in-cone size using the same rule as the rate panel.
    score, best_shape = quops_score_from_summary(plot_df.assign(passes=plot_df["decision"].eq("Pass")))
    if best_shape is not None:
        score_label = "Gated QUOPS score" if "score_passes" in plot_df else "Largest certified tested size" if "multiplicity" in plot_df and plot_df["multiplicity"].eq("hochberg").all() else "Largest tested pass"
        ax.scatter([best_shape[1]], [best_shape[0]], s=165, marker="D", facecolors="#4B7600", edgecolors="#203400", linewidths=1.3, label=f"{score_label}: {score:,} QUOPS ({best_shape[0]} qubits)", zorder=6).set_gid("quops-score")
    else:
        ax.text(0.03, 0.97, "Gated score: 0\nNo gated score pass" if "score_passes" in plot_df else "Tested score: 0\nNo passing in-cone shape", transform=ax.transAxes, ha="left", va="top", fontsize=10, color="#555555").set_gid("quops-score-label")
    ax.legend(loc="lower right", framealpha=0.95)
    return fig


def _quops_scan_options(summary_df, threshold):
    """Validate a scan's table format and explicitly declared test settings."""
    hochberg_results = "multiplicity" in summary_df and summary_df["multiplicity"].eq("hochberg").any()
    if hochberg_results and not summary_df["multiplicity"].eq("hochberg").all():
        raise ValueError("a scan cannot mix Hochberg and other test families")
    gated_results = "multiplicity" in summary_df and summary_df["multiplicity"].eq("gated_hochberg").any()
    if gated_results and not summary_df["multiplicity"].eq("gated_hochberg").all():
        raise ValueError("a scan cannot mix gated and other test families")
    required = {"width", "size", "decision"}
    if not required.issubset(summary_df.columns):
        raise ValueError(f"summary_df is missing required columns: {sorted(required - set(summary_df.columns))}")
    threshold_values = pd.to_numeric(summary_df["threshold"], errors="coerce").dropna().unique() if "threshold" in summary_df.columns else np.array([])
    if threshold is None and len(threshold_values) > 1:
        raise ValueError("a scan plot requires a common threshold")
    threshold = _validate_threshold(float(threshold_values[0]) if threshold is None and len(threshold_values) == 1 else QUOPS_THRESHOLD if threshold is None else threshold)
    if len(threshold_values) and not np.allclose(threshold_values, threshold, rtol=0.0, atol=1.0e-12):
        raise ValueError("plot threshold must match the threshold used for shape decisions")
    test_description = "Individual 95% tests"
    if hochberg_results:
        alpha_values = pd.to_numeric(summary_df["alpha"], errors="raise").unique()
        if len(alpha_values) != 1:
            raise ValueError("a Hochberg plot requires one common family error budget")
        test_description = f"Hochberg tests; family error budget = {_probability('alpha', alpha_values[0]):.0%}"
    elif gated_results:
        required = {"score_passes", "score_alpha", "region_alpha", "score_confidence", "region_confidence"}
        if not required.issubset(summary_df.columns):
            raise ValueError(f"gated results are missing columns: {sorted(required - set(summary_df.columns))}")
        budgets = {}
        for column in ("score_alpha", "region_alpha"):
            values = pd.to_numeric(summary_df[column], errors="raise").unique()
            if len(values) != 1:
                raise ValueError("a gated plot requires common score and region error budgets")
            budgets[column] = _probability(column, values[0])
        if sum(budgets.values()) >= 1.0:
            raise ValueError("score and region error budgets must sum to less than one")
        for column, expected in (("score_confidence", 1.0 - budgets["score_alpha"]), ("region_confidence", 1.0 - sum(budgets.values()))):
            if not np.allclose(pd.to_numeric(summary_df[column], errors="raise"), expected, rtol=0.0, atol=1.0e-12):
                raise ValueError("gated confidence metadata does not match the separate error budgets")
        test_description = f"Gated score: {1.0 - budgets['score_alpha']:.0%}; combined region: {1.0 - sum(budgets.values()):.0%} confidence"
    if hochberg_results or gated_results:
        if "inference_method" in summary_df:
            methods = summary_df["inference_method"].dropna().unique()
            if len(methods) != 1:
                raise ValueError("a notebook scan requires one sampling-inference method")
            if methods[0] not in ("hoeffding", "normal"):
                raise ValueError("unknown sampling-inference method")
            test_description += "; finite-sample component-mean bounds" if methods[0] == "hoeffding" else "; normal approximation"
    return threshold, test_description


def plot_quops_scan(summary_df: pd.DataFrame, rate_df: pd.DataFrame | None = None, *, title: str = "QUOPS capability", label: str | None = None, pass_column: str = "passes", threshold: float | None = None, shape_figsize: tuple[float, float] = (9.0, 4.8), dpi: int = 140, ax=None, assume_monotonic: bool = False):
    """Plot the notebook's shape decisions and optional simulator rates.

    Accept width, size and decision columns, including a gated Not tested tail.
    Hochberg and gated results retain their declared threshold and budgets.
    The utility cone is geometric. Set assume_monotonic=True to draw the
    capability region's downward closure (SI Eq. 68). This assumes mean
    polarization cannot increase with width or size; tests do not verify it.
    Pass ax to draw a single capability panel in an existing figure.

    Supply rate_df from calculate_quops_rate to return one figure with two
    panels. Rates match decisions by (width, size). The diamonds select the
    largest passing in-cone size and its rate. Gated tables select the score
    from score_passes alone, even if larger region shapes pass.
    Omit ax for this layout.
    """
    import matplotlib.pyplot as plt
    threshold, test_description = _quops_scan_options(summary_df, threshold)
    if not isinstance(assume_monotonic, (bool, np.bool_)):
        raise ValueError("assume_monotonic must be a boolean")
    if assume_monotonic and ("multiplicity" not in summary_df or not summary_df["multiplicity"].isin(["hochberg", "gated_hochberg"]).all()):
        raise ValueError("inferred capability regions require family-controlled results")

    # Match the separately calculated rates to the tested shapes.
    if rate_df is not None:
        if ax is not None:
            raise ValueError("omit ax when supplying rate_df; the combined scan creates two axes")
        rate_columns = {"width", "size", "omega_simulated"}
        if not rate_columns.issubset(rate_df.columns):
            raise ValueError(f"rate_df is missing required columns: {sorted(rate_columns - set(rate_df.columns))}")
        rate_summary = summary_df.drop(columns=["omega_simulated"], errors="ignore").merge(rate_df[["width", "size", "omega_simulated"]], on=["width", "size"], how="left", validate="one_to_one")
        rate_summary["decision"] = rate_summary["decision"].astype(str).str.lower()
        rate_summary[pass_column] = rate_summary[pass_column].fillna(False).eq(True) & rate_summary["decision"].eq("pass") if pass_column in rate_summary else rate_summary["decision"].eq("pass")
        fig, axes = plt.subplots(1, 2, figsize=(2 * shape_figsize[0], shape_figsize[1]), dpi=dpi, layout="constrained")
        _plot_quops_individual_scan(summary_df, title=title, label=label, threshold=threshold, figsize=shape_figsize, dpi=dpi, ax=axes[0], test_description=test_description, assume_monotonic=assume_monotonic)
        rate_values = pd.to_numeric(rate_summary["omega_simulated"], errors="coerce")
        if (np.isfinite(rate_values) & rate_values.gt(0.0)).any():
            plot_quops_rate(rate_summary, pass_column=pass_column, label=label, dpi=dpi, ax=axes[1])
        else:
            axes[1].set_title("Simulator QUOPS throughput")
            axes[1].set_xlabel("Circuit size, s (QUOPS)")
            axes[1].set_ylabel("Simulator throughput (QUOPS/s)")
            axes[1].text(0.5, 0.5, "No positive finite rates available", transform=axes[1].transAxes, ha="center", va="center")
        return fig
    return _plot_quops_individual_scan(summary_df, title=title, label=label, threshold=threshold, figsize=shape_figsize, dpi=dpi, ax=ax, test_description=test_description, assume_monotonic=assume_monotonic)


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
    eligible = source_df["decision"].astype(str).str.casefold().eq("pass") if "decision" in source_df else pd.Series(True, index=source_df.index)
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
            rate_text = f"{selected_rate:,.2f}" if selected_rate >= 1.0 else f"{selected_rate:.3g}"
            _highlight_quops_point(ax, best_shape[1], selected_rate, f"Rate at score\n{rate_text} QUOPS/s", gid="quops-rate-at-score", legend_label=f"Rate at ${score_symbol}={score}$")
        else:
            ax.text(0.97, 0.04, f"Rate at ${score_symbol}={score}$ unavailable", transform=ax.transAxes, ha="right", fontsize=9, color="#555555")
    else:
        ax.text(0.97, 0.04, "Rate at score unavailable\nNo gated score pass" if "score_passes" in source_df else "Rate at score unavailable\nNo passing in-cone shape", transform=ax.transAxes, ha="right", fontsize=9, color="#555555")
    ax.set_xlabel("Circuit size, $s$ (QUOPS)", fontsize=11)
    heading = "Simulator QUOPS throughput" if rate_column == "omega_simulated" else "QUOPS rate by shape"
    ax.set_ylabel("Simulator throughput (QUOPS/s)" if rate_column == "omega_simulated" else "Rate (QUOPS/s)", fontsize=11)
    ax.set_title(heading if label is None else f"{heading}\n{label}", fontsize=11)
    if rate_column == "omega_simulated":
        ax.text(0.97, 0.96, "BR/RR sampling-call proxy\nREF excluded", transform=ax.transAxes, ha="right", va="top", fontsize=9, color="#555555")
    if omitted:
        ax.text(0.01, 0.01, f"{omitted} invalid/nonpositive rate row(s) omitted", transform=ax.transAxes, fontsize=8, color="#555555", va="bottom")
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.9, fontsize=9)
    if standalone:
        fig.tight_layout()
    return fig


# ===========================================================================
# Acquisition and provenance: seeds, raw counts and resumable archives
# ===========================================================================

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


def save_notebook_evidence(data: pd.DataFrame, measurements: Iterable[dict], shapes: Iterable[dict], settings: dict, output=None) -> Path:
    """Save a new notebook evidence directory after checking raw-count consistency.

    Each measurement has width, size, circuit_index, mirror_index, and
    counts/targets dictionaries keyed by br, rr and ref. A source circuit
    manifest is optional; when supplied, its consistency is checked.
    The directory contains the observed table, raw records, exact local source
    snapshots and a settings/environment manifest. Existing paths are rejected.
    This records evidence; it does not establish MCFE model validity.
    """
    shapes = default_shape_schedule(shapes)
    validate_acquisition(data, shapes, num_circuits=settings["num_circuits"], num_mirrors=settings["num_mirrors"], shots=settings["shots"])
    keys = ["width", "size", "circuit_index", "mirror_index"]
    observed = {tuple(row[key] for key in keys): row for row in data.to_dict("records")}
    records = list(measurements)
    seen = set()
    circuits = {}
    for record in records:
        key = tuple(_nonnegative_integer(name, record[name]) for name in keys)
        if key not in observed or key in seen:
            raise ValueError("raw measurement has an unexpected or duplicate circuit/mirror identity")
        seen.add(key)
        row = observed[key]
        if "circuit" in record:
            circuit = record["circuit"]
            length = len(circuit["gates"])
            if circuit["width"] != row["width"] or any(len(circuit[name]) != length for name in ("q0s", "q1s", "angles", "phis", "lams")) or circuit_size(circuit["gates"]) != row["size"]:
                raise ValueError("saved source circuit does not match its observed width/size or gate arrays")
            circuit_json = canonical_json(circuit)
            if key[:3] in circuits and circuits[key[:3]] != circuit_json:
                raise ValueError("mirrors of one original circuit have different saved source circuits")
            circuits[key[:3]] = circuit_json
        for ensemble in ("br", "rr", "ref"):
            histogram = counts_dict(record["counts"][ensemble])
            polarization = effective_polarization_from_counts(histogram, row["width"], record["targets"][ensemble])
            if sum(histogram.values()) != row[f"shots_{ensemble}"]:
                raise ValueError(f"saved {ensemble} counts do not reproduce their observed shot count")
            if not math.isclose(polarization, float(row[f"lambda_{ensemble}"]), rel_tol=1.0e-12, abs_tol=1.0e-12):
                raise ValueError(f"saved {ensemble} counts/target do not reproduce their observed effective polarization")
    if seen != set(observed):
        raise ValueError("raw measurements do not cover every observed circuit/mirror")

    # Validate all records before creating an exclusive directory; no prior evidence is overwritten.
    created = datetime.now(timezone.utc)
    destination = Path(output) if output is not None else Path("results") / f"notebook_{created.strftime('%Y%m%dT%H%M%S')}_{uuid4().hex[:8]}"
    destination.mkdir(parents=True, exist_ok=False)
    source_directory = destination / "sources"
    source_directory.mkdir()
    sources = {}
    for source in (Path(__file__), Path(__file__).with_name("quops.ipynb")):
        if source.exists():
            contents = source.read_bytes()
            (source_directory / source.name).write_bytes(contents)
            sources[source.name] = sha256(contents).hexdigest()
    data.to_csv(destination / "observations.csv", index=False)
    with (destination / "measurements.jsonl").open("x") as stream:
        for record in records:
            stream.write(canonical_json(record) + "\n")
    configuration = {"schema_version": 1, "depth_unit": "paper", "shapes": shapes, "settings": json_safe(settings), "software": software_manifest(), "source_sha256": sources, "timing_scope": "sum_of_M1_M2_sample_call_spans", "is_hardware_rate": False}
    manifest = {"created_at_utc": created.isoformat(), "configuration": configuration, "configuration_sha256": sha256(canonical_json(configuration).encode()).hexdigest(), "files_sha256": {name: file_sha256(destination / name) for name in ("observations.csv", "measurements.jsonl")}, "observations": len(data), "raw_counts_verified": True}
    atomic_write_json(destination / "manifest.json", manifest)
    return destination


def validate_runtime() -> dict:
    """Check the optimization-boundary feature before collecting any data."""
    try:
        kernel = build_kernel([], [], [], [], 1)
    except TypeError as exc:
        raise RuntimeError('This CUDA-Q build lacks atomic quantum regions. Use requirements.txt and the README setup instructions; do not silently disable MCFE boundaries.') from exc
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
    metrics = {"Resource requests": counts.operation_counts.get("resource_request", 0), "Syndrome extractions": counts.syndrome_rounds, "Peak live blocks (including injection)": counts.patches_peak, "Carriers in peak block templates": counts.patches_peak * code.block.size}
    return _show(metrics, table, title=f"P2: {code.name} encoded resources", x="instruction", y="count", xlabel="Static Fabric instruction instances (log scale)", color="#76B900")


# ===========================================================================
# Public notebook imports
# ===========================================================================

__all__ = ['BR_kernel', 'QUOPS_THRESHOLD', 'REF_kernel', 'RR_kernel', 'RunArchive', 'adjoint_C', 'atomic_write_json', 'calculate_gamma_hat_ws', 'calculate_gamma_lower', 'calculate_mean_lambdas', 'calculate_p_value', 'calculate_quops_rate', 'calculate_sigma', 'canonical_json', 'counts_dict', 'default_shape_schedule', 'describe_noise', 'effective_polarization_from_counts', 'enumerate_shapes', 'file_sha256', 'gated_hochberg_test', 'hamming_distance', 'hochberg_test', 'in_utility_cone', 'json_safe', 'kernel_manifest', 'make_depolarizing_noise', 'make_helios1_noise', 'normalized_mcfe_polarization', 'omega_simulated', 'plot_p0_comparison', 'plot_p1_capacity', 'plot_qec_overhead', 'plot_quops_rate', 'plot_quops_scan', 'plot_quops_timing', 'quops_rate', 'quops_score_from_summary', 'randomized_compilation', 'save_notebook_evidence', 'sample_mcfe_caps', 'sample_quops_arrays', 'sample_quops_circuit', 'shape_size', 'show_p0', 'show_p2', 'software_manifest', 'utility_cone_bounds', 'validate_acquisition', 'validate_runtime']
