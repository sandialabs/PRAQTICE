# What is my quantum computer good for? Quantum capability learning with physics-aware neural networks

This repository is the temporary supplemental material of [What is my quantum computer good for? Quantum capability learning with physics-aware neural networks](https://arxiv.org/abs/.). 

The repository contains copies of the code, data sets, models, and notebooks used to create the simulation results in the paper. To reproduce the experimental results, you will need to access the data [here](https://zenodo.org/records/7829489) and edit the notebooks provided herein accordingly. The directory is organized as follows:
   - This README.md file.
   - The environment.yml file used to set up the conda environment used throughout the paper.
   - The `ml` folder contains the bulk of the code. Inside are the files:
      - `tools.py` and `newtools.py` which contain many helper functions for processing the data,
      - `newneuralnets.py` which contains the quantum-physics-aware network classes,
      - `modeltools.py` which contains several functions used when instantiating the NNs,
      - and `encoding.py` which contains many of the functions used to encode a `pygsti` circuit into a tensor.
   - The `paper-simulations-fidelity` subdirectory. Inside are:
      - Notebooks used to generate the simulation data, train the qpa-NNs and CNNs, tune the CNNs, and generate plots.
      - Subdirectories `experiment_{experiment_id}` that contain the data used in each simulation trial.
   - A subdirectory containing additional plots not found in the paper.

## Requirements

To install requirements:

```setup
conda env create -f environment.yml
```

You will also need to manually install [pyGSTi](https://github.com/sandialabs/pyGSTi) version 0.9.11.2.post173. Instructions are available [here](https://github.com/sandialabs/pyGSTi).

## Datasets

To generate a dataset that is similar to one of the 10 simulated datasets used in the paper, follow the instructions in the following notebooks located in the simulations-fidelity subdirectory:

```dataset generation
1-circuit-generation.ipynb
2-simulations.ipynb
3-package-circuits.ipynb
```

The first notebook (1-circuit-generation.ipynb) will generate a collection of random *i.i.d.*-layer circuits or random mirror circuits. The second notebook (2-simulations.ipynb) will generate (or load) a pygsti error model, and simulate the collection of random circuits to compute their entanglement fidelities. The third notebook (3-package-circuits.ipynb) will process each pair of circuits and fidelities into the tensor format required for model training. The third notebook will also partition the dataset into training, validation, and testing subsets. For simplicity, circuits are packaged as I(c) + P(c) + S(c).

## Training

To train the model(s) used in the simulation section of the paper, follow the instructions in the following notebooks (also located in the simulations-fidelity subdirectory):

```train
4-training-fidelity.ipynb
5-training-cnn.ipynb
```

The first notebook (4-training-fidelity.ipynb) creates a new instance of the `CircuitErrVec` class (the `keras` implementation of the physics-aware network used in the simulation section of the paper). The first notebook also trains the model, saves the model weights, and its predictions, either on the random *i.i.d.*-layer circuits or the random mirror circuits. The second notebook (5-training-cnn.ipynb) does the same, except using a hyperparameter-tuned CNN.

## Evaluation

To evaluate a model use either of the following two notebooks:

```eval
6-paper-plots.ipynb
7-data-analysis.ipynb
```

The first notebook (6-paper-plots.ipynb) is a copy of the notebook used to generate the plots in the paper. Within the first notebook are instructions for generating: (i) scatter plots of the physics-aware and convolutional networks' predictions, (ii) scatter plots of the physics-aware networks' prediction errors; and (iii) summary statistics (MAE) and histograms for every model's predictions. 

## Pre-trained Models

You can find the weights of each pre-trained physics-aware neural network, and the pre-trained hyperparameter-tuned CNNs in the models subdirectory within each dataset's directory (i.e., experiment_{experiment_id}).

## Results

Our models achieve the following performance across the 10 simulated datasets used in the paper:

| Model name         |  Circuit type  | Mean absolute error  |    Std. dev.    |       LLR       |
| ------------------ | -------------  | -------------------- | --------------- | --------------- |
|                    | random  i.i.d. |         .18%         |       ***       |       ***       |
| Physics-aware NN   | -------------- | -------------------- | --------------- | --------------- |
|                    |     mirror     |         .72%         |       ***       |       ***       |
| ------------------ | -------------- | -------------------- | --------------- | --------------- |
|                    | random  i.i.d. |         .75%         |       ***       |       ***       |
| CNN                | -------------- | -------------------- | --------------- | --------------- |
|                    |     mirror     |         2.6%         |       ***       |       ***       |
| ------------------ | -------------- | -------------------- | --------------- | --------------- |
 

Below is the summary figure (Fig. 2) from the paper's experimental section.

![A scatter plot of prediction errors (a) on the `ibmq_vigo` dataset and (b) a plot depicting the (i) physics-aware networks', (ii) state-of-the-art CNNs'; and (iii) fine-tuned CNNs' prediction accuracy on each experimental dataset](plots/combined-experimental-plot.png)


## License and copyright

This code is meant for release in `pygsti`, which is copyrighted and operates under an Apache 2.0 license. Please DO NOT distribute without permission from the author.
