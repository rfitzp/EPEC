# NEOCLASSICAL

## Description

   Program to read FLUX data from fFile(s), profile data from pFile(s) 
   and cFiles(s), and write neoclassical data needed by PHASE to nFile
   - Stage 3:
     Read input files, calculate neoclassical data at rational 
     surfaces, and output nFile

## Requirements

   - blitz++ library (https://github.com/blitzpp/blitz)
   - gsl library (https://www.gnu.org/software/gsl)
   - BLAS library (http://www.netlib.org/blas)
   - NETCDF library (http://unidata.ucar.edu)

## Contents

### /Documentation
- Neoclassical.tex: Latex description of program algorithm
	  
### /Source

- Makefile: GNU makefile
- *.f90: Fortran_90 source files
- *.h: C++ header files
- *.cpp: C++ source files
	 
### /Inputs

- Neoclassical.nml: Fortran_90 namelist control file
- fFile: Data from PHASE
- /fFiles: Contains fFiles for interpolation
  - Index: List of fFile names and experimental times
  - fFiles: Actual fFiles
- pFile: Profile data
- /pFiles: Contains pFiles for interpolation
  - Index: List of pFile names and experimental times
  - pFiles: Actual pFiles
- cFile: Perpendicular diffusivity data
- /cFiles: Contains cFiles for interpolation
  - Index: List of cFile names and experimental times
  - cFiles: Actual cFiles
	  
### /Outputs

- nFile: File containing data for PHASE
- /nFiles: Directory containing interpolated nFiles 
  - Index: List of nFile names and experimental times
  - nFiles: Actual nFiles
- Stage3.nc: NETCDF data file from Stage 3 calculation
	  
### /Plots

- /Stage3
  -*.py: Python scripts to plot Stage 3 data

### pFile Format

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
	
	n "psinorm nb(10^20/m^3) dpbdpsiN"
	for (int i = 0; i < n; i++)
	PSI, NB, dNBdPSI

	n "psinorm pb(kPa) dnbdpsiN"
	for (int i = 0; i < n; i++)
	PSI, PB, dPBdPSI

	n "psinorm ptot(kPa) dptotdpsiN"
	for (int i = 0; i < n; i++)
	PSI, PTOT, dPTOTdPSI
	
	n "psinorm omeg(kRad/s) domeg/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WPHI, dWPHIdPSI

	n "psinorm omegp(kRad/s) domegp/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WTHE, dWTHEdPSI

	n "psinorm omegvb(kRad/s) domevb/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WVB, dWVBdPSI

	n "psinorm omegpp(kRad/s) domepp/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WPP, dWPPdPSI
	
	n "psinorm omgeb(kRad/s) domgeb/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WEB, dWEBdPSI

	n "psinorm er(kV/m) der/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, ER, dERdPSI

	n "psinorm ommvb(kRad/s) dommvb/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WMVB, dWMVBdPSI

	n "psinorm ommpp(kRad/s) dommpp/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WMPP, dWMPPdPSI

	n "psinorm omevb(kRad/s) domevb/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WEVB, dWEVBdPSI

	n "psinorm omepp(kRad/s) domepp/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WEPP, dWEPPdPSI

	n "psinorm kpol(km/s/T) dkpol/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, KPOL, dKPOLdPSI

	n "psinorm omghb() domghb/dpsiN"
	for (int i = 0; i < n; i++)
	PSI, WMGB, dWMGBdPSI
	
	n "psinorm nz1(10^20/m^3) dnz1dpsiN"
	for (int i = 0; i < n; i++)
	PSI, NI, dNIdPSI

	n "psinorm vtor1(km/s) dvtor1psiN"
	for (int i = 0; i < n; i++)
	PSI, VTOR1, dVTOR1dPSI

	n "psinorm vpol1(km/s) dvpol1psiN"
	for (int i = 0; i < n; i++)
	PSI, VPOL1, dVPOL1dPSI
	
	n "N Z A of ION SPECIES"
	for (int i = 0; i < n; i++)
	N, Z, A (i=0 impurity, i=1 majority; i=2 fast)

 - PSI:   Normalized poloidal flux
 - NE:    Electron number density (10^20/m^3)
 - TE:    Electron temperature (keV)
 - NI:    Thermal ion number density (10^20/m^3)
 - TI:    Thermal ion temperature (keV)
 - NB:    Fast ion number density (10^20/m^3)
 - WPHI:  Impurity ion toroidal angular velocity on outboard midplane (krad/s)
 - WTHE:  Impurity ion toroidal angular velocity on outboard midplane (krad/s)	
 - WEB:   ExB frequency (krad/s)
 - NI :   Impurity ion number density (10^20/m^3)
 - N:     Ion atomic number
 - Z:     Ion charge (units of e)
 - A:     Ion mass number

 *Fields can occur in any order. Additional fields are ignored.*
 
 ### cFile Format
 
    n
    for (int i = 0; i < n; i++)
    PSI, CHI_PHI, CHI_E, D_PERP, CHI_I
	
 - PSI:     Normalized poloidal flux	
 - CHI_PSI: Perpendicular toroidal momentum diffusivity (m^2/s)
 - CHI_E:   Perpendiocular electron energy diffusivity (m^2/s)
 - D_PERP:  Perpendicular particle diffusivity (m^2/s)
 - CHI_I:   Perpendiocular ion energy diffusivity (m^2/s)
