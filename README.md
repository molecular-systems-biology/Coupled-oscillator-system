# Coupled-oscillator-system
This repository contains the code necessary to estimate the parameters for Kuramoto model of cell cycle and metabolic oscillators. Parameter estimation part is written in GAMS and the simulation part is written in MATLAB. 
- `input.gdx`: This is the measured input to the optimization problem which consists for phase differences and frequencies. It is in GDX format which is the standard data input format of GAMS. 
- `parameter_estimation.gms`: This is the GAMS script to perform the multistart optimizations to estimate the parameters of the coupled oscillators system. 
- `Kuramoto_Simulations.mlx`: It is a MATLAB Live Script file that contains the code steady-state simulations with the estimated parameter sets. 
- `parray.mat`: MAT-file that contains the estimated parameter sets, phase differences and frequencies. 
- `kuramoto.m`: Contains the MATLAB function for Kuramoto model, it is called from `Kuramoto_Simulations.mlx` to integrate the model equations.
