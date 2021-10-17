WAVEFORM
========

Keywords
--------

namelist

Short Description
-----------------

Description of WAVEFORM namelist

Waveform Control
----------------

TYPE
  Waveform type (1=programmed, 2=spike, 3=repeated ramp)


Type 1 Waveform
---------------

NCTRL 
  Number of RMP control points in waveform
TCTRL 
  List of RMP control times (ms)
ICTRL 
  List of RMP currents at control times (kA)
PCTRL 
  List of relative phases of RMP coil currents at control times (units of PI)

Type 2 Waveform
---------------

SSTART
  Time of start of RMP spike (ms)
SEND 
  Time of end of RMP spike (ms)
SAMP
  RMP coil currents at start of RMP spike (kA)
SPHA
  Relative phase of RMP coil currents at start of spike (units of PI)
SPVE
  Phase velocity of RMP coil currents during spike (krad/s)
BACK
  Background RMP coil currents (kA)

Type 3 Waveform
---------------

RPERIOD
  Repeated RMP ramp repetition period (ms) 
RSTART 
  RMP coil currents at start of ramp (first ramp starts at t=0) (kA)
REND 
  RMP coil currents at end of ramp (kA)
RPHA 
  Relative phase of RMP coil currents during ramp (units of PI)

NTM Trigger
-----------

MPOL
 Poloidal mode number of NTM whose triggering is under investigation 
AMIN
 Minimum amplitude of applied type 2 RMP spike (kA)
AMAX
 Maximum amplitude of applied type 2 RMP spike (kA)
PSTART
 Minumum duration of applied type 2 RMP spike (starts at SSTART) (ms)
PEND 
 Maximum duration of applied type 2 RMP spike (ms)
DP 
 Difference between durations of applied type 2 RMP spikes (ms) 
