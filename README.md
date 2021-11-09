# Generation

## How to Make a Circuit Model

The routines for building SBOL circuit components are in `builders.py`
Import the module and run via a script like `make_sbol_models.py`

## Generating LaTeX Equations

The routines for generating LaTeX from SBOL circuits are in `latex_generation.py`
Import the module and run via a script like `sbol_to_latex.py`

## Generating MATLAB Code

The routines for generating LaTeX from SBOL circuits are in `matlab_generation.py`
Import the module and run via a script like `sbol_to_matlab.py`


# Curve fitting

The function `parameter_fitting.m` is a Matlab function to fit an ODE model of a
Cre-ON kill switch from to the data on said circuit from Chylinski et al.
It uses the ODE in `fit_Cre_on_Kill_Switch.m`.
