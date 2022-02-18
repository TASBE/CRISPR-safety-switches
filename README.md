# CRISPR Time-Delayed Safety Switch Models

[![gh-action badge](https://github.com/TASBE/CRISPR-safety-switches/actions/workflows/python-app.yml/badge.svg)](https://github.com/TASBE/CRISPR-safety-switches/actions)

This repository contains supplementary information for the paper: 

>  Helen Scott, Dashan Sun, Jacob Beal, and Samira Kiani. 
  "Simulation-Based Engineering of Time-Delayed Safety Switches for Safer Gene Therapies", 
  _to appear_.

If you make use of any of the contents of this repository, please cite this paper.

## Abstract
CRISPR-based gene editing is a powerful tool with great potential for applications in the treatment of many inherited and acquired diseases. 
The longer that CRISPR gene therapy is maintained within a patient, however, the higher the likelihood that it will result in problematic side effects such as off-target editing or immune response. 
One approach to mitigating these issues is to link the operation of the therapeutic system to a safety switch that autonomously disables its operation and removes the delivered therapeutics after some amount of time. 
We present here a simulation-based analysis of the potential for regulating the time-delay of such a safety switch using one or two transcriptional regulators and/or recombinases. 
Combinatorial circuit generation identifies 30 potential architectures for such circuits, which we evaluate in simulation with respect to tunability, sensitivity to parameter values, and sensitivity to cell-to-cell variation. 
This modeling predicts one of these circuit architectures to have the desired dynamics and robustness, which can be further tested and applied in the context of CRISPR therapeutics.

## Contents of Repository

This repository contains:
- Python code that generates SBOL models each potential safety switch architecture
- Python code to export SBOL models into ODE equations in LaTeX and executable ODE simulation models in Matlab.
- Regression tests for validating the operation of the model generators.
- Simulation runners for experimentation with the models.
- Results and figures generated from these simulations.


## Running the Code

### How to Make a Circuit Model

The routines for building SBOL circuit components are in `builders.py`
Import the module and run via a script like `make_sbol_models.py`

### Generating LaTeX Equations

The routines for generating LaTeX from SBOL circuits are in `latex_generation.py`
Import the module and run via a script like `sbol_to_latex.py`

### Generating MATLAB Code

The routines for generating LaTeX from SBOL circuits are in `matlab_generation.py`
Import the module and run via a script like `sbol_to_matlab.py`

## Curve fitting

The function `parameter_fitting.m` is a Matlab function to fit an ODE model of a
Cre-ON kill switch from to the data on said circuit from Chylinski et al.
It uses the ODE in `fit_Cre_on_Kill_Switch.m`.
