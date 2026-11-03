# QUOPS on quantum hardware

- Preserve the benchmark definitions: source size counts rotations plus twice the CNOT count; use threshold `1/√e` and score shapes within `w² ≤ size ≤ w³`. [§§II–III](https://arxiv.org/html/2609.12146v1#S2)
- Use simulations or pilot runs to locate the pass boundary; test sizes more densely near it. [§VIII.A](https://arxiv.org/html/2609.12146v1#S8.SS1)
- Optimize layout, routing and native compilation before constructing randomized references. Preserve the original QUOPS circuit as the fidelity target. [§II.C](https://arxiv.org/html/2609.12146v1#S2.SS3)
- Preserve compilation barriers between mirror halves and around preparation/readout caps. [§VII.D.1](https://arxiv.org/html/2609.12146v1#S7.SS4)
- Preserve reference-circuit layers through final compilation, including internal layers. Justify any architecture-required exceptions. [§II.C](https://arxiv.org/html/2609.12146v1#S2.SS3)
- Use fresh preparation, reference and measurement randomizations; prefer per-shot updates where the backend supports them. [§VII.D.1](https://arxiv.org/html/2609.12146v1#S7.SS4)
- Validate bit order, mirror targets, preserved barriers and the native reference compiler's MCFE assumptions before the full run.
- Keep REF observations independent with the current APIs. Shared references require different uncertainty handling.
- Separate capability measurement from throughput optimization. Increase shots per circuit to reduce switching overhead, then confirm the selected shape still passes. [§VIII.B, D](https://arxiv.org/html/2609.12146v1#S8)
- Measure throughput at the score-selected shape using first-shot-to-last-shot QPU timestamps. Include required control, waveform compilation, reset and calibration overhead; exclude unrelated jobs. State whether REF time is included. [§VII.A.4](https://arxiv.org/html/2609.12146v1#S7.SS1), [§VIII.D](https://arxiv.org/html/2609.12146v1#S8.SS4)
- If postselecting, use a validated acceptance rule and retain total and accepted shot counts. Report ordinary and postselected results separately and include the discarded-shot cost in the rate. [Eq. (72)](https://arxiv.org/html/2609.12146v1#S4.SS5)
- Save circuits, raw counts, targets, seeds, compiler settings, calibration records and device timestamps.
