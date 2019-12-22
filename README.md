REMD_analysis
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/REMD_analysis.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/REMD_analysis)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/REMD_analysis/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/REMD_analysis/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/REMD_analysis/branch/master)

## Description
`REMD_analysis` is a Python package of data analysis tools for replica exchange molecular dynamics (REMD) simulations.

## Installation and Testing
All the Python scripts in this package are written in Python 3. Currently the package can be installed by following the commands below:
```
git clone https://github.com/wehs7661/REMD_analysis.git
cd REMD_analysis
pip install -e .
```
To perform the unit tests and functional tests (still underdevelopment) of this package, run:
```
python test_REMD_analysis.py    ; unit tests
python test_REMD_analysis.sh    ; functional tests
```

## Usage

#### 1. Prepartion of simulation input files 
This package provides a simple way to prepare simulations input files (such as `.mdp` and `.tpr` files for each replica) for Hamiltonian replica exchange simulations and submitting jobs to Summit or Bridges. By running the command `bash HREMD.sh`, the following promopts will be invoked for the user to specify relevant paramters:
```
This shell script prepares all the files needed to run a Hamiltonian replica exchange, including .mdp files and .tpr 
files, given the common .gro file and a template .mdp file. A job submission script will also be generated.
Please input the name of the job: 
Please input the prefix of the file names: 
Please input the number of replicas: 
Please input the number of requested nodes: 
Please input the simulation time (in hour): 
Will the job be submitted to Summit or Bridges?
```
The script for preparing temperature replica exchange simulations are currently under development.

#### 2. `REMD_analysis.py`: Visualization of transition/overlap matrices and state-time plot 


## External Links
To get more realized about the theory and the implmentation of replica exchange simulation, we recommend the following materials:

- Tutorials
- Literature (provided in the folder `REMD_analysis/papers`)
  - Temperature replica exchange
  - Hamiltonian replica exchange
  - Optimization of replica exchange simulations
  - Multistate Bennett Acceptance Ratio (MBAR)


## Copyright

Copyright (c) 2019, Wei-Tse


## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
