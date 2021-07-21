# FLUX

## Description

   Program to read gFile(s) and write equilibrium data needed by NEOCLASSICAL and PHASE to fFile.
   - Stage 1:
     Read gFile and output equilibrium data for Stage 2
   - Stage 2:
     Construct flux coordinate system. Calculate metric quantities. Locate rational surfaces.
     Caculate neoclassical data. Calculate tearing stability matrix. 
	 
## Requirements

   - blitz++ library (https://github.com/blitzpp/blitz)
   - gsl library (https://www.gnu.org/software/gsl)
   - BLAS library (http://www.netlib.org/blas)
   - asymptote (for plots) (https://asymptote.sourceforge.io)
   
## Contents

### /Documentation

- Flux.tex: Latex description of program algorithm
	  
### /Source

- Makefile: GNU makefile
- *.f90: Fortran_90 source files
- *.h: C++ header files
- *.cpp: C++ source files
- /Test: Test routine for toroidal functions
	 
### /Inputs

- Flux.nml: Fortran_90 namelist control file
- gFile: Equilibrium gFile
- /gFiles: Contains equilibrium gFiles for interpolation
  - Index: List of gFile names and experimental times
  - gFiles: Actual gFiles
  - rFile: rFile
	  
### /Outputs

- fFile: File containing data for NEOCLASSICAL and PHASE
- /fFiles: Directory containing fFiles associated with interpolated gFiles
 - Index: List of fFile names and experimental times
 - fFiles: Actual fFiles
- /Stage1: Data files from Stage 1 calculation
- /Stage2: Data files from Stage 2 calculation
	  
### /Plots

- /Stage1
 - *.asy: Asymptote scripts to plot Stage 1 data
- /Stage2
  - *.asy: Asymptote scripts to plot Stage 2 data

## gFile Format

    read (100, '(a48, 3i4)') string, i3, NRBOX, NZBOX
    read (100, '(5e16.9)') RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero
    read (100, '(5e16.9)') RAXIS, ZAXIS, PSIAXIS, PSIBOUND, B0
    read (100, '(5e16.9)') CURRENT, zero, zero, zero, zero
    read (100, '(5e16.9)') zero, zero, zero, zero, zero
 
    read (100, '(5e16.9)') (T (i), i = 1, NRBOX)
    read (100, '(5e16.9)') (P (i), i = 1, NRBOX)
    read (100, '(5e16.9)') (TTp (i), i = 1, NRBOX)
    read (100, '(5e16.9)') (Pp (i), i = 1, NRBOX)
  
    read (100, '(5e16.9)') ((PSI (i, j), i = 1, NRBOX), j = 1, NZBOX)

    read (100, '(5e16.9)') (Q (i), i = 1, NRBOX)

    read (100, '(2i5)') NBOUND, NLIM

    read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
    read (100, '(5e16.9)') (RLIM (i), ZLIM (i), i = 1, NLIM)
  
  *Anything after this is ignored*
  
## rFile Format

    RPLUS RMINUS ZPLUS ZMINUS
	
  *First line is lFile, second is uFile, third is mFile*	
