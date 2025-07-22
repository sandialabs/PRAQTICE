# Efficient Two Qubit GST

In this subrepository you can find files associated with the analysis of the performance of single and two qubit GST under a number of different
experiment designs. This includes different choices of germ selection techniques, and different methods for performing fiducial pair reduction.
This analysis is related to the paper "Two-Qubit Gate Set Tomography with Fewer Resources" which can be found here: https://arxiv.org/pdf/2307.15767

The general structure of this subrepository is as follows:

- data: In this directory can be found all of the data sets and serialized output necessary to perform the performance analysis of the one- and two-qubit experiment designs considered in https://arxiv.org/pdf/2307.15767. This includes directories for each of the serialized experiment designs, a sub-directory of pandas dataframes containing relevant analysis results, and a number of compressed numpy arrays with pre-computed fisher information matrices for each experiment design. This data can be found in a separate repository located TBD and should be downloaded and copied into this directory before running any code which relies on it.
- numerical_analysis : In this directory all of the scripts and shell scripts used to perform the performance analysis of each of the experiment designs can be found.
- paper_figures : This directory contains jupyter notebooks specifically relevant to the generation of the experiment design visualization methods highlighted and included in the aforementioned paper.

Additionally, two pip requirements files are provided with specific tested configurations. Note that there are two files as the one- and two-qubit analyses contained within this subrepository were originally performed in two slightly different configuration environments, and as such may not be 100% compatible with each other.