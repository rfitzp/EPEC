PHASE
=====

Short Description
-----------------

Program to perform EPEC island dynamics simulation.

Keywords
--------

island, RMP

Long Description
-----------------
 
Program to read FLUX data from fFile(s), NEOCLASSICAL data from 
nFile(s), GPEC data from uFile(s), mFiles(s), and lFile(s), 
and perform island dynamics simulation. 
Final state of plasma saved in sFile.

Stage 4:
  Read data and calculate vacuum island widths as function of relative 
  RMP coil phase
Stage 5:
  Perform island dynamics simulation in fixed equilibrium
Stage 6:
  Perform island dynamics simulation in dynamic equilibrium  

Contents
--------
	  
Inputs:
 Phase.nml:
  Fortran_90 namelist control file
 Waveform.nml:
  Fortran_90 namelist file specifing RMP coil current waveform
 fFile:
  Equilibrium data from FLUX
 fFiles:
  Directory containing fFiles for interpolation
    Index:
     List of fFile names and experimental times
    fFiles:
     Actual fFiles 
 nFile:
  Neoclassical/profile data from NEOCLASSICAL
 nFiles:
  Directory containing nFiles for interpolation
    Index:
     List of nFile names and experimental times
    nFiles:
     Actual nFiles
 lFile:
   Ideal response data to 1kA current flowing in lower RMP coil set from GPEC
 lFiles:
   Directory containing lFiles for interpolation
    Index:
     List of lFile names and experimental times
    lFiles:
     Actual lFiles
 uFile:
   Ideal response data to 1kA current flowing in (optional) upper RMP coil set from GPEC
 uFiles:
   Directory containing uFiles for interpolation
    Index:
     List of uFile names and experimental times
    uFiles:
     Actual uFiles
 mFile:
   Ideal response data to 1kA current flowing in (optional) middle RMP coil set from GPEC
 mFiles:
   Directory containing mFiles for interpolation
    Index:
     List of mFile names and experimental times
    mFiles:
     Actual mFiles
	  
Outputs:
 sFile:
  File specifying final plasma state (used for restarting calculation)
 Stage4.nc:
  NETCDF file from Stage 4 calculation
 Stage5.nc:
  NETCDF file from Stage 5 calculation
 ncFiles/\*.nc:
  NETCDF files from Stage 6 calculation  
	  
Plots:
  \*.py:
   Python scripts to plot Stage 4 and 5 data
   
  ../EPEC/PLOTS/\*.py:
   Python scripts to plot Stage 6 data 
  
  
lFile, uFile, mFile Format
--------------------------

Files are *gpec_singfld_n?.out* files from GPEC with singular field and ascii flags set
