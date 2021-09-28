

EPEC
====

Short Description
-----------------
Suite of programs to simulate island dynamics in presence of RMP

Keywords
--------
island, dynamics, RMP, GPEC

Long Description
----------------
Suite of programs to simulate multi-harmonic magnetic island chain dynamics 
in presence of resonant magnetic perturbation in time-varying toroidal 
tokamak equilibrium.

Relevant Publications
---------------------
* R. Fitzpatrick, SangKyeun Kim, and Jaehyun Lee.
  *Modeling of q95 windows for the suppression of edge localized modes by resonant magnetic perturbations in the KSTAR tokamak*.
  Phys. Plasmas **28**, 082511 (2021).
* R. Fitzpatrick.
  *Further modeling of q95 windows for the suppression of edge localized modes by resonant magnetic perturbations in the DIII-D tokamak*.
  Phys. Plasmas **28**, 022503 (2021).
* R. Fitzpatrick.
  *Modeling q95 windows for the suppression of edge localized modes by resonant magnetic perturbations in the DIII-D tokamak*.
  Phys. Plasmas **27**, 102511 (2020).
* R. Fitzpatrick, and A.O. Nelson.
  *An improved theory of the response of DIII-D H-mode discharges to static resonant magnetic perturbations and its implications for the suppression of edge localized modes*.
  Phys. Plasmas **27**, 072501 (2020).

Weblinks
--------
* `EPEC Homepage <http://farside.ph.utexas.edu/img/EPEC-documentation/epec.html>`_
* `Asymptotic Matching Approach to Modeling Tearing Mode Dynamics in Tokamak Plasmas <http://farside.ph.utexas.edu/talks/ifs2021.pdf>`_
* `Predicting Stability Windows in q95 for RMP-Induced ELM Suppression in H-Mode Tokamak Discharges <http://farside.ph.utexas.edu/talks/PoPWebinar.pdf>`_

Contents
--------
/FLUX 
  Program to read gFile(s) and write equilibrium data needed by 
  NEOCLASSICAL and PHASE to fFile(s)

/NEOCLASSICAL 
  Program to read FLUX data from fFile(s), profile data from pFile(s) 
  and cFiles(s), and write neoclassical data needed by PHASE to nFile(s)

/PHASE 
  Program to read FLUX data from fFile(s), NEOCLASSICAL data from 
  nFile(s), GPEC data from uFile(s), mFiles(s), and lFile(s), and 
  perform island dynamics simulation. Final 
  state of plasma saved in sFile.

/RESCALE 
  Program to rescale equilibrium gFile (and pFile) so as to modify q95 by modifying 
  toroidal plasma current while keeping toroidal magnetic field the same

/GPEC
  Independent OMFIT module used to calculate ideal-MHD response of tokamak
  plasma equilibirum to resonant magnetic perturbation
