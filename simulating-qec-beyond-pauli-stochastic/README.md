This directory contains code to reproduce the simulations in the paper "Simulating Quantum Error Correction beyond Pauli Stochastic Errors", https://arxiv.org/abs/2603.18457. This paper develops a method for generating a DEM from a sparse error generator model for physical-qubit gates and uses it to study the performance of FTQC primitives. The code for the DEM generation method is in pyGSTi, and this repo contains code specific to the simulation case studies. 

Subdirectories:
cultivation/ Code for simulating an analyzing Cliffordized magic state cultivation circuits (constructing "sensitivity matrices"). cutlivation_sweep_h_error.py runs the DEM construction and Monte Carlo simulations, h_sweep_analysis.ipynb plots the results, and sensitivty_analysis_for_cultivation.ipynb produces sensitivity matrices for cultivation circuits. 

memory_demos/
    gross_code/ Code for simulating memory circuits for the [[144,12,12]] code and comparing arbitrary CPTP error models with their Pauli twirled analogues. Contains two scripts for running simulations (one for parallelized simulations), a notebook for producing plots (gross_code_simulation_plots.ipynb), and a notebook demonstrating the simulation workflow (bb_code_error_prop_workflow_demo.ipynb)
    
surface_code/ Code for simulating surface code memory with (a) random CPTP error models (surface_random_models.py), and (b) S+H error models with fixed CNOT generator infidelity (surface_coherent_error_vs_stochastic_error.ipynb), and an accomanying plotting notebook (paper_threshold_plots.ipynb).

validation/ Code for reproducing the validation tests for small circuits using exact simulation methods from LoQS and Qiskit.
    loqs_simulation_analysis.ipynb: Analysis comparing DEMs to weak simulation of the true model with LoQS
    steane_code_dem_baseline.ipynb: Code determining the DEM events that occur with a standard circuit-level depolarizing noise model
    h_error_simulation_data.ipynb: Plotting comparisons between our DEMs and state vector simulation for H-only error models
    test_random_h_models_sweep_p.py: Script to simulate surface code syndrome extraction with random crosstalk-free coherent error models
    random_h_models_pauli_twirl.py: Script to estimate detection history probabilities for Pauli twirled versions of the random H-only error models
    pauli_twirled_tools.py: Helper functions for analysis of Pauli twirled models. 
    loqs_dems_helpers.py: Tools for interacting with the LoQS simulations
    loqs_test_scha.py: Script to generate dems for Steane code syndrome extraction with random sparse CPTP error models, and simulate those models with LoQS

The code here requires a decoder for the Gross code. Here, we use the Beam Search decoder, available at https://github.com/ionq-publications/BeamSearchDecoder (commit 084a475b05fb64308103317a1ce5a0c4b0be58aa). Similar results can be obtained by replacing this decoder with other BP-based decoders, such as BP-OSD. 

Additionallly, we use helper functions from Gidney et al.'s magic state cultivation code, available at https://zenodo.org/records/13777072.