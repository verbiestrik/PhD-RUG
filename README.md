# Adaptive simulation of moment models 

This repository contains the code I developed for the model-adaptive simulation of moment models with application in free-surface flow modelling and rarefied gas dynamics.

**Overview**

The software was developed as part of my PhD project in numerical and computational mathematics at the University of Groningen. A common workflow is first deriving moment model PDEs, then solving the moment model PDEs numerically, and finally comparing the results with those of a reference solution. 

**Features**

- Automatic derivation of two-dimensional shallow water moment models and one-dimensional Grad-type moment models for rarefied gases
- Numerical simulation of moment model PDEs that have a hyperbolic structure
- Numerical simulation of reference equations for one-dimensional shallow flow and radially symmetric flow, adapted from the code developed in https://github.com/ShallowFlowMoments/Supplements2018
- Numerical simulation of linear kinetic moment models with stability analysis

**Installation**

git clone https://github.com/verbiestrik/PhD-RUG.git

cd PhD-RUG 

pip install -r requirements.txt

**Usage**

For numerically solving non-linear moment model PDEs, navigate to SpatiallAdaptiveMomentModels/Nonlinear-systems. Next, choose the desired inputs in the Config-files/config.txt file. Then run python main.py.

**Requirements**

Python 3.11+
Dependencies listed in requirements.txt

Add your name or a way to reach you, if relevant.
