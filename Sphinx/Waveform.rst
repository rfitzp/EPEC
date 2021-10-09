WAVEFORM
========

Keywords
--------

namelist

Short Description
-----------------

Description of WAVEFORM namelist

Namelist Variables
------------------

TYPE
  Waveform type (1=programmed, 2=spike, 3=repeated ramp) 
NCTRL 
  Number of RMP control points in type 1 programmed waveform
TCTRL 
  List of RMP control times (type 1) (ms)
ICTRL 
  List of RMP currents at control times (type 1) (kA)
PCTRL 
  List of relative phases of RMP coil currents at control times (type 1) (units of PI)
SSTART
  Time of start of RMP spike (type 2) (ms)
SEND 
  Time of end of RMP spike (type 2) (ms)
WAMOD
  Amplitude modulation frequency of RMP coil currents during RMP spike (type 2) (krad/s)
WPMOD
  Relative phase velocity of RMP coil currents during RMP spike (type 2) (krad/s)
SAMP
  RMP coil currents at start of RMP spike (type 2) (kA)
SPHA
  Relative phase of RMP coil currents at start of spike (type 2) (units of PI)
BACK
  Background RMP coil currents (type 2) (kA)
RPERIOD
  Repeated RMP ramp repetition period (type 3) (ms) 
RSTART 
  RMP coil currents at start of ramp (first ramp starts at t=0) (type 3) (kA)
REND 
  RMP coil currents at end of ramp (type 3) (kA)
RPHA 
  Relative phase of RMP coil currents during ramp (type 3) (units of PI)
  
