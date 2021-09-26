NEOCLASSICAL NAMELIST
=====================

Short Description
-----------------

Description of NEOCLASSICAL namelist

Keywords
--------

namelist

Namelist Variables
------------------

EXB
  Switch for calculating ExB frequency from neoclassical theory:
   0:
    Use ExB frequency profile from pFile
   1:
    Use ExB frequency profile derived from pFile toroidal and poloidal velocity profiles (should be same as EXB = 0)
   2:
    Use ExB frequency profile derived from pFile toroidal velocity profile and neoclassical poloidal velocity profile
IMPURITY 
  Flag for inclusion of impurities in calculation
NEUTRAL 
  Flag for inclusion of neutrals in claculation
NTYPE 
  If 0/1 then flux-surface averaged neutral density distribution exponential/Lorentzian
NN 
  Flux-surface averaged neutral density on LCFS (m^-3)
LN 
  Flux-surface averaged neutral density decay lengthscale (m)
YN 
  Neutral peaking factor at X-point
SVN
  Charge exchange rate constant (m^3/s)
EN 
  Incoming neutral energy/ion energy
INTF
  Flag for fFile interpolation
INTP
  Flag for pFile interpolation
INTC
  Flag for cFile interpolation
CATS 
  Flag for linear-only interpolation of cFiles
DMIN
  Minimum allowed value of diffusvities at resonant surfaces (m^2/s)
DMAX
  Maximum allowed value of diffusvities at resonant surfaces (m^2/s)  
TIME 
  Experimental time (ms) (only relevant with interpolation of pFiles/fFiles/cFiles)
TAUMIN
  Minimum allowed value of tau (parameter in calculation of linear layer widths)
COULOMB
  Coulomb logarithm
NSMOOTH
  Number of smoothing cycles for calculation of higher derivatives of profiles

