FLUX NAMELIST
=============

Short Description
-----------------

Description of FLUX namelist

Keywords
--------

namelist

Namelist Variables
------------------

NTOR 
  Toroidal mode number of RMP
MMIN
  Minimum poloidal mode number of resonant surfaces included in calculation
MMAX 
  Maximum poloidal mode number of resonant surfaces included in calculation
PSILIM
  Maximum PsiN for safety-factor calculation
PSIRAT
  Resonant surfaces in region PsiN > PSIRAT are ignored
PSIPED
  PsiN value at top of pedestal
INTG 
  Flag for interpolating gFiles
TIME 
  Experimental time (ms) (only relevant when INTG = 1)
RW
  Radius of resistive wall (units of minor radius)
NPSI 
  Number of points in PsiN grid
PACK
  Packing index for PsiN grid:
   PsiN[i] = PSILIM * (1. - (1. - s)^PACK) for i = 0, NPSI-1, where s = i /(NPSI-1)
NTHETA
  Number of points in theta grid (must be odd)
NNC
  Number of neoclassical harmonics included in calculation
NSMOOTH 
  Number of smoothing cycles for calculation of higher derivatives of q
H0 
  Initial integration step-length
ACC 
  Adaptive integration accuracy
ETA
  Regularization factor for Green's function
NEOANG 
  Flag for using neoclassical, as opposed to straight, angle in E-matrix calculation
