# QUOPS with CUDA-Q

This directory provides CUDA-Q examples for **QUOPS (quantum universal operation performance system)**, described in [*Benchmarking the computational power of quantum computers*](https://arxiv.org/html/2609.12146v1) by Timothy Proctor et al. (arXiv:2609.12146v1, 2026). The notebook implements the paper's mirror-circuit fidelity estimation (MCFE) workflow, estimates passing circuit sizes and simulator throughput, and explores logical implementation costs through P0 synthesis, P1 placement and P2 Steane encoding.

| File or directory | Purpose |
| --- | --- |
| [quops.ipynb](quops.ipynb) | Main tutorial: BR/RR/REF measurements, polarization and Hochberg tests, capability/rate plots, and logical resource examples. |
| [src.py](src.py) | Circuit generation, randomized compilation, inference, plotting and optional acquisition/calibration command-line tools. |
| [requirements.txt](requirements.txt) | Pinned direct dependencies and Jupyter tools; pip installs their dependencies automatically. |
| [plots/](plots/) | Saved physical capability/rate, synthesis-precision and code-memory figures displayed and explained in the notebook. |

## Run the notebook

Use **Linux x86-64 with Python 3.12.3** to match the recorded environment. The default notebook uses CPU simulation and does not require GPU execution. From this directory:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip==26.2.1
python -m pip install -r requirements.txt
python -m pip check
python -m ipykernel install --sys-prefix --name quops --display-name "Python (QUOPS)"
python -m jupyterlab quops.ipynb
```

Select **Python (QUOPS)** and run the notebook from the first cell. Keep `src.py` and `plots/` beside it. For command-line usage, run these in the same environment from this directory:

```bash
python src.py --help
python src.py calibrate --help
```

## Environment and saved results

`requirements.txt` pins the direct dependencies and Jupyter tools to the versions used in the **12 September 2026** environment that successfully executed all nine notebook code cells. It uses `cudaq==0.16.0`, which installs the matching CUDA-Q Logical package on Linux and selects the CUDA binary distribution automatically. Transitive dependencies can resolve to newer compatible versions; these pins are not a complete environment lock.

Running the small notebook example generates new measurements. The saved figures illustrate a physical scan up to 30 qubits, synthesis precision and code-memory resources. Their raw measurements, original environment records and acquisition scripts are not included in this directory, so running the demo does not reproduce those larger studies. Recorded simulator rates depend on the machine and timing conditions; the logical resource examples do not measure a logical QUOPS score.
