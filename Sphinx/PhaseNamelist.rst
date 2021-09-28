PHASE NAMELIST
==============

Short Description
-----------------

Description of PHASE namelist

Keywords
--------

namelist

Namelist Variables
------------------

PMAX
  RMP coil relative phase scan from -PMAX*M_PI/2 to PMAX*M_PI/2 in Stage 4
STAGE5
  Flag for proceeding to Stage5 calculation
TSTART
  Simulation start time (ms)
TEND 
  Simulation end time (ms)
DT 
  Data recorded every DT ms
MID
  Number of coil sets:
   1:
    Requires lFile
   2:
    Requires lFile and uFile
   3:
    Requires lFile, uFile, and mFile
COPT
  Switch for optimaiztion of relative amplitudes/phases of coils currents:
   0:
    No optimization
   1:
    Maximize drive at resonant surface closest to top of pedestal (restricted)
   2:
    Maximize drive at resonant surface closest to top of pedestal (unrestricted)
   3:
    Maximize drive at resonant surface closest to top of pedestal and minimize drive at innermost resonant surface (unrestricted)
CORE
  Core drive minimization factor when COPT = 3 (0.0 = no minimization, 1.0 = complete minmization)
FREQ
  Switch for selecting natural frequency type:
   0:
    Use linear/nonlinear natural frequency with linear layer width as switch
   1:
    Use linear/nonlinear natural frequency with electron pressure flattening width as switch 
   2:
    Use linear/ExB/nonlinear natural frequency with linear layer width as switch
   3:
    Use w_natural = FFAC * w_linear + (1-FFAC) * w_EB
FFAC
  Parameter for selecting natural frequency (when FREQ = 3)
LIN 
  Flag for purely linear island dynamics simulation
CXD 
  Flag for including charge exhange damping in angular equations of motion
BOOT
  Flag for including bootstrap current in Rutherford equations
CURV
  Flag for including magnetic field-line curvature in Rutherford equations
POLZ
  Flag for including ion polarization current in Rutherford equations
TAUW
  Time constant of resistive wall (ms)
CHIR
  Maximum Chirikov parameter for vacuum islands
INTF
  Flag for fFile interpolation
INTN
  Flag for nFile interpolation
INTU
  Flag for uFile/mFile/lFile interpolation
NATS
  Flag for linear-only nFile interpolation
RATS
  Flag for linear-only interpolation of uFiles/mFiles/lFiles
OLD 
  Flag for restarting old calculation
HIGH
  Flag for higher-order transport analysis
SCALE
  GPEC scalefactor
NFLOW
  Number of flow harmonics included in calculation

