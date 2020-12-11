# ISLANDDYNAMICS

## Description

   Program to simulate multi-harmonic magnetic island dynamics in presence of static, externally generated, 
   resonant magnetic  perturbation in time-evolving toroidal tokamak discharge by calling
   the FLUX, NEOCLASSICAL, and PHASE codes.
   - Stage 6: Perform simulation and output data.
	
## Requirements

   - blitz++ library (https://github.com/blitzpp/blitz)
   - gsl library (https://www.gnu.org/software/gsl)
   - BLAS library (http://www.netlib.org/blas)
   - asymptote (for plots) (https://asymptote.sourceforge.io)
   
## Contents

### /Source

- Makefile: GNU makefile
- *.f90: Fortran_90 source files
- *.cpp: C++ source files
	 
### /Inputs

- Island.in: Fortran_90 namelist file
- Waveform.in: Fortran_90 namelist file specifying RMP coil current waveform
	  
### /Outputs

 - /Stage6: Data files from Stage 6 calculation
	  
### /Plots

- /Stage6
  - *.asy: Asymptote scripts to plot Stage 6 data
