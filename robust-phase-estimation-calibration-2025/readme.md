This is the readme file for the supplemental material to be used in conjunction with the
paper "Heisenberg-limited calibration of entangling gates with robust phase estimation",
which can be found here:  https://arxiv.org/abs/2502.06698.

    Table of Contents
    1.  Folder contents
    2.  Installation instructions

1.  Folder contents

This folder contains the following files:

    README.txt : This file.
    
    requirements.txt : Requirements file listing necessary python packages.
    
    CZ_3_phase_RPE_demo.ipynb : Jupyter notebook demonstrating how to use robust phase
    estimation for estimating the three rotation angles in a two-qubit controlled-phase 
    gate.

2.  Installation instructions

We recommend creating a fresh virtual environment as to not potentially interfere with
any environments you may already be using.  Instructions here use venv but other options
(e.g., conda) should work just fine.  Terminal instructions are as follows:

cd /path/to/robust-phase-estimation-calibration-2025

python -m venv rpe

source rpe_venv/bin/activate

pip install --upgrade pip

pip install requirements.txt

pip install ipykernel

python -m ipykernel install --user --name=rpe

Note that this will install both pyRPE, which performs the RPE analysis, and pygsti which 
is used for quantum circuit construction and dataset manipulation.  Only the basic 
installation of pygsti is included.  For more advanced pygsti usage (e.g., two-qubit gate 
set tomography, consider a complete installation of pygsti; see pygsti.info.)
