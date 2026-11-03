# QUOPS with CUDA-Q

This directory implements the mirror-circuit fidelity estimation (MCFE) workflow in [*Benchmarking the computational power of quantum computers*](https://arxiv.org/html/2609.12146v1). The main notebook estimates QUOPS circuit sizes and simulator throughput.

| File or directory | Purpose |
| --- | --- |
| [quops.ipynb](quops.ipynb) | Main workflow: BR/RR/REF measurements, polarization, shape tests and plots. |
| [src.py](src.py) | Circuit generation, randomized compilation, analysis, plotting and result storage. |
| [requirements.txt](requirements.txt) | Pinned direct dependencies and Jupyter tools. |
| [HARDWARE_GUIDE.md](HARDWARE_GUIDE.md) | Best practices for adapting and running QUOPS on quantum hardware. |
| [logical_quops.ipynb](logical_quops.ipynb) | Optional synthesis, resource estimation and ideal Steane execution example. |
| [plots/](plots/) | Saved example figures. |

## Run the notebook

Use **Linux x86-64 with Python 3.12**. The default target is `density-matrix-cpu`; no GPU is required. From this directory:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip check
python -m ipykernel install --sys-prefix --name quops --display-name "Python (QUOPS)"
python -m jupyterlab quops.ipynb
```

Select **Python (QUOPS)** and run the cells in order. Keep `src.py` and `plots/` beside the notebook. Let acquisition finish before running analysis. Raw counts, targets, observations and settings are saved under `results/`.


## Logical example

`logical_quops.ipynb` defines a CUDA-Q `BR` kernel that calls the forward circuit, randomized inverse and caps. It imports that kernel directly with `cql.compiler.import_cudaq`, then lowers it through P0 synthesis, P1 placement, P2 QEC realization and P3 physical compilation. 

