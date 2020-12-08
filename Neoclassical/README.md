# NEOCLASSICAL

## Description

   Program to read FLUX data from fFile(s), profile data from pFile(s) and cFiles(s), and write neoclassical data needed by PHASE to nFile
   - Stage 3:
     Read input files, calculate neoclassical data at rational surfaces, and output nFile
   
## Contents

 ### /Documentation
    - Neoclassical.tex
      Latex description of program algorithm
	  
 ### /Source
     - Makefile
	 GNU Makefile
	- *.f90
	 Fortran_90 source files
	- *.h
	 C++ header files
	- *.cpp
	 C++ source files
	 
 ### /Inputs
	- Neoclassical.in
	  Fortran_90 namelist control file
	- fFile
	  Data from PHASE
	- pFile
	  Profile data
	- cFile
	  Perpendicular diffusivity data
	  
 ### /Outputs
    - nFile
	  File containing data for PHASE
	- /nFiles
      Directory containing precalculated nFiles
	- /Stage3
	  Data files from Stage 3 calculation
	  
 ### /Plots
    - /Stage3
	 -- README
	    Description of Asymptote scripts
	 -- *.asy
	    Asymtptote scripts to plot Stage 3 data

## pFile Format

	n "psinorm ne(10^20/m^3) dnedpsiN"
	for (int i = 0; i < n; i++)
	PSI, NE, dNEdPSI
	
	n "psinorm te(KeV) dtedpsiN"
	for (int i = 0; i < n; i++)
	PSI, TE, dTEdPSI
	
	n "psinorm ni(10^20/m^3) dnidpsiN"
	for (int i = 0; i < n; i++)
	PSI, NI, dNIdPSI
	
	n "psinorm ti(KeV) dtidpsiN"
	for (int i = 0; i < n; i++)
	PSI, TI, dTIdPSI
	
	n "psinorm nb(10^20/m^3) dnbdpsiN"
	for (int i = 0; i < n; i++)
	PSI, NB, dNBdPSI
	
	n "psinorm omeg(kRad/s) domeg/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WPHI, dWPHIdPSI
	
	n "psinorm omgeb(kRad/s) domgeb/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WEB, dWEBdPSI
	
	n "psinorm nz1(10^20/m^3) dnz1dpsiN"
	for (int i = 0; i < n; i++)
	PSI, NI, dNIdPSI
	
	n "N Z A of ION SPECIES"
	for (int i = 0; i < n; i++)
	N, Z, A
	
 *Fields can occur in any order. Additional fields are ignored.*
 
 ## cFile Format
 
    n
    for (int i = 0; i < n; i++)
    PSI, CHI_PHI, CHI_E, D_PERP
