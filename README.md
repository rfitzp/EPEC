# EPEC
Program to simulate multi-harmonic magnetic island dynamics in presence of resonant magnetic perturbation 
in time-varying toroidal tokamak equilibrium.

RESCALE        - rescales equilibrium gFile to modify q_95, keeping Btor the same

FLUX           - reads gFile(s) and generates equilibrium data needed by NEOCLASSICAL and PHASE

NEOCLASSICAL   - reads pFile(s) and generates neoclassical data needed by PHASE

PHASE          - reads GPEC data (uFile(s) and lFile(s)) and performs island dynamics simulation in fixed equilibrium

ISLANDDYNAMICS - performs island dynamics simulation in time-varying equilibrium

