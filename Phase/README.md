# PHASE

## Description
 
 Program to read FLUX data from fFile(s), NEOCLASSICAL data from 
 nFile(s), GPEC data from uFile(s), mFiles(s), and lFile(s), 
 and perform island dynamics simulation in fixed equilibrium. 
 Final state of plasma saved in sFile.
 - Stage 4:
	  Read data and calculate vacuum island widths as 
	  function of RMP coil phase
 - Stage 5:
	  Perform island dynamics simulation

## Requirements

   - blitz++ library (https://github.com/blitzpp/blitz)
   - gsl library (https://www.gnu.org/software/gsl)
   - BLAS library (http://www.netlib.org/blas)
   - NETCDF library (http://unidata.ucar.edu)
   
## Contents

### /Documentation

- Phase.pdf: Description of program algorithm
	  
### /Source

- Makefile: GNU makefile
- *.f90: Fortran_90 source files
- *.h: C++ header files
- *.cpp: C++ source files
	 
### /Inputs

- Phase.nml: Fortran_90 namelist control file
- Waveform.nml: Fortran_90 namelist file specifing RMP 
  coil current waveform
- fFile: Data from FLUX
- /fFiles: Directory containing fFiles for interpolation
  - Index: List of fFile names and experimental times
  - fFiles: Actual fFiles
- nFile: Data from NEOCLASSICAL
- /nFiles: Directory containing nFiles for interpolation
  - Index: List of nFile names and experimental times
  - nFiles: Actual nFiles
- uFile: GPEC data for upper RMP coil set
- /uFiles: Directory containing uFiles for interpolation
  - Index: List of uFile names and experimental times
  - uFiles: Actual uFiles
- mFile: GPEC data for (optional) middle RMP coil set
- /mFiles: Directory containing mFiles for interpolation
  - Index: List of mFile names and experimental times
  - mFiles: Actual mFiles
- uFile: GPEC data for lower RMP coil set
- /uFiles: Directory containing uFiles for interpolation
  - Index: List of uFile names and experimental times
  - mFiles: Actual uFiles
	  
### /Outputs

- sFile: File specifying final plasma state (used for 
  restarting calculation)
- Stage4.nc: NETCDF data files from Stage 4 calculation
- Stage5.nc: NETCDF data files from Stage 5 calculation
	  
### /Plots

- *.py: Python scripts to plot Stage 4 and Stage 5 data		

## uFile, mFile, lFile Format

 Files are gpec_singfld_n?.out files from GPEC with
 singular field and ascii flags set
