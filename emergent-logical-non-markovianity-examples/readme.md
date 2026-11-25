This is the readme file for the supplemental material to be used in conjunction with the
paper "Emergent Non-Markovianity in Logical Qubit Dynamics",
which can be found here:  https://arxiv.org/abs/2512.08893

    Table of Contents
    1.  Folder contents
    2.  Installation instructions

1.  Folder contents

This folder contains the following files:

    readme.md : This file.
    
    requirements.txt : Requirements file listing necessary python packages.
    
    EmergentLogicalNonMarkovianityExamples.ipynb : Jupyter notebook which provides examples 
    of emergent logical non-Markovianity in the 3-qubit repetition ([3,1,3]) code and the 5-qubit perfect ([[5,1,3]]) code.

2.  Installation instructions

We recommend creating a fresh virtual environment as to not potentially interfere with
any environments you may already be using.  Instructions here use conda but other options
(e.g., venv) should work just fine.  Terminal instructions are as follows:

cd /path/to/emergent-logical-non-markovianity-examples

conda create --name ELNME

conda install pip

pip install requirements.txt
