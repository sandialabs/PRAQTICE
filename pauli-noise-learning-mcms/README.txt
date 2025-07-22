Code from the paper "Pauli noise learning for mid-circuit measurements"
by Jordan Hines and Timothy Proctor 
https://arxiv.org/abs/2406.09299

This work introduces mid-circuit measurement cycle benchmarking (MCM-CB), a benchmarking method that estimates the fidelity of a randomly compiled mid-circuit measurement. In the above paper, we introduce a theory for this method, and for Pauli noise learning of MCMs more broadly, and we demonstrate MCM-CB in simulations with up to 10 qubits and on IBM Q processors with up to 4 qubits. We show how to use the data from this method to estimate the probabilities of individual Pauli error rates, and do a detailed analysis of errors during a single-qubit MCM on an IBM Q device. This folder includes the code to reproduce our results. 

Subdirectories:
ibmq: Contains notebooks to run MCM-CB on an IBM Q processor. NOTE: Requires an IBM Q account and device access to run. Produces Fig. 1b-c and Fig. 3 in the main text. 

simulations: Contains scripts to run the simulations of MCM-CB shown in the main text and SM and notebooks to reproduce all plots. 
Each notebook 