RESCALE
=======

Short Description
-----------------

Program to rescale equilibrium and profile data needed by EPEC.

Keywords
--------

rescale, equilibrium

Long Description
-----------------

Program to read gFile, pFile, and cFile, rescale equilibrium,
and write out new gFile, pFile, and cFile
	 
Contents
--------

/Inputs:
  Rescale.nml:
   Fortran_90 namelist control file
  gFile:
   Initial gFile
  pFile:
   Initial pFile
  cFile:
   Initial cFile
	  
/Outputs:
  gFile:
   Rescaled gFile
  pFile:
   Rescaled pFile
  cFile
   Rescaled cFile
	  
gFile Format
------------

| read (100, '(a48, 3i4)') string, i3, NRBOX, NZBOX
| read (100, '(5e16.9)') RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF
| read (100, '(5e16.9)') RAXIS, ZAXIS, PSIAXIS, PSIBOUND, B0
| read (100, '(5e16.9)') CURRENT, zero, zero, zero, zero
| read (100, '(5e16.9)') zero, zero, zero, zero, zero
 
| read (100, '(5e16.9)') (T (i), i = 1, NRBOX)
| read (100, '(5e16.9)') (P (i), i = 1, NRBOX)
| read (100, '(5e16.9)') (TTp (i), i = 1, NRBOX)
| read (100, '(5e16.9)') (Pp (i), i = 1, NRBOX)
  
| read (100, '(5e16.9)') ((PSI (i, j), i = 1, NRBOX), j = 1, NZBOX)
 
| read (100, '(5e16.9)') (Q (i), i = 1, NRBOX)

| read (100, '(2i5)') NBOUND, NLIM
 
| read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
| read (100, '(5e16.9)') (RLIM (i), ZLIM (i), i = 1, NLIM)
 
*Anything after this is ignored*
  
pFile Format
------------

| n "psinorm ne(10^20/m^3) dnedpsiN"
| for (int i = 0; i < n; i++)
| PSI, NE, dNEdPSI
	
| n "psinorm te(KeV) dtedpsiN"
| for (int i = 0; i < n; i++)
| PSI, TE, dTEdPSI
	
| n "psinorm ni(10^20/m^3) dnidpsiN"
| for (int i = 0; i < n; i++)
| PSI, NI, dNIdPSI
	
| n "psinorm ti(KeV) dtidpsiN"
| for (int i = 0; i < n; i++)
| PSI, TI, dTIdPSI
 
| n "psinorm nb(10^20/m^3) dpbdpsiN"
| for (int i = 0; i < n; i++)
| PSI, NB, dNBdPSI

| n "psinorm pb(kPa) dnbdpsiN"
| for (int i = 0; i < n; i++)
| PSI, PB, dPBdPSI

| n "psinorm ptot(kPa) dptotdpsiN"
| for (int i = 0; i < n; i++)
| PSI, PTOT, dPTOTdPSI
	
| n "psinorm omeg(kRad/s) domeg/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WPHI, dWPHIdPSI

| n "psinorm omegp(kRad/s) domegp/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WTHE, dWTHEdPSI

| n "psinorm omegvb(kRad/s) domevb/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WVB, dWVBdPSI

| n "psinorm omegpp(kRad/s) domepp/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WPP, dWPPdPSI
	
| n "psinorm omgeb(kRad/s) domgeb/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WEB, dWEBdPSI

| n "psinorm er(kV/m) der/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, ER, dERdPSI

| n "psinorm ommvb(kRad/s) dommvb/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WMVB, dWMVBdPSI

| n "psinorm ommpp(kRad/s) dommpp/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WMPP, dWMPPdPSI

| n "psinorm omevb(kRad/s) domevb/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WEVB, dWEVBdPSI

| n "psinorm omepp(kRad/s) domepp/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WEPP, dWEPPdPSI

| n "psinorm kpol(km/s/T) dkpol/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, KPOL, dKPOLdPSI

| n "psinorm omghb() domghb/dpsiN"
| for (int i = 0; i < n; i++)
| PSI, WMGB, dWMGBdPSI
	
| n "psinorm nz1(10^20/m^3) dnz1dpsiN"
| for (int i = 0; i < n; i++)
| PSI, NI, dNIdPSI

| n "psinorm vtor1(km/s) dvtor1psiN"
| for (int i = 0; i < n; i++)
| PSI, VTOR1, dVTOR1dPSI

| n "psinorm vpol1(km/s) dvpol1psiN"
| for (int i = 0; i < n; i++)
| PSI, VPOL1, dVPOL1dPSI
	
| n "N Z A of ION SPECIES"
| for (int i = 0; i < n; i++)
| N, Z, A (i=0 impurity, i=1 majority; i=2 fast)

PSI:
  Normalized poloidal flux
NE:
  Electron number density (10^20/m^3)
TE:
  Electron temperature (keV)
NI:
  Thermal ion number density (10^20/m^3)
TI:
  Thermal ion temperature (keV)
NB:
  Fast ion number density (10^20/m^3)
WPHI:
  Impurity ion toroidal angular velocity on outboard midplane (krad/s)
WTHE:
  Impurity ion toroidal angular velocity on outboard midplane (krad/s)	
WEB:
  ExB frequency (krad/s)
NI:
  Impurity ion number density (10^20/m^3)
N:
  Ion atomic number
Z:
  Ion charge (units of e)
A:
  Ion mass number

*Fields can occur in any order. Additional fields are ignored.*

cFile Format
------------
 
| n
| for (int i = 0; i < n; i++)
| PSI, CHI_PHI, CHI_E, D_PERP, CHI_I
	
PSI:
 Normalized poloidal flux	
CHI_PSI:
 Perpendicular toroidal momentum diffusivity (m^2/s)
CHI_E:
 Perpendicular electron energy diffusivity (m^2/s)
D_PERP:
 Perpendicular particle diffusivity (m^2/s)
CHI_I:
 Perpendicular ion energy diffusivity (m^2/s)
