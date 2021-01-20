# WINDOWFIND

## Description

   Program to find RMP suppression windows in q95 versus Irmp space
   by calling ISLANDDYNAMICS code.
	 
## Requirements

   - blitz++ library (https://github.com/blitzpp/blitz)
   - gsl library (https://www.gnu.org/software/gsl)
   - BLAS library (http://www.netlib.org/blas)
   - asymptote (for plots) (https://asymptote.sourceforge.io)
   
## Contents
	  
### /Source

  - Makefile: GNU makefile
  - *.f90: Fortran_90 source files
  - *.h: C++ header files
  - *.cpp: C++ source files
	 
### /Inputs

  - Window.nml: Fortran_90 namelist control file
  
### /Outputs

 - *.txt: Output files
	  
### /Plots

 - *.asy: Asymptote scripts to plot RMP ELM suppresion windows
