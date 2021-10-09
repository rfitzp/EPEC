EPEC NAMELIST
=============

Short Description
-----------------

Description of EPEC namelist

Keywords
--------

namelist

Namelist Variables
------------------

TSTART
  Simulation start time (ms)
TEND
  Simulation end time (ms)
DTR
  Time interval between successive RESCALE calculations (ms)
DTF
  Time interval between successive fFile/nFile calculations (ms)  
DTP
  Time interval between successive PHASE runs (ms)
RESTART
   Flag for restarting each PHASE run, rather than using previous run as initial condition
ANSTART
  Initial AN for RESCALE n_e scan
ANEND
  Final AN for RESCALE n_e scan
ATSTART
  Initial AT for RESCALE T_e scan
ATEND
  Final AT for RESCALE T_e scan
ARSTART
  Initial AR for RESCALE R_0 scan
AREND
  Final AR for RESCALE R_0 scan
ACSTART
  Initial AC for RESCALE CHI_PHI scan
ACEND
  Final AC for RESCALE CHI_PHI scan       
APSTART
  Initial AP for RESCALE P scan
APEND
  Final AP for RESCALE P scan      
AWSTART
  Initial AW for RESCALE w_E scan
AWEND
  Final AW for RESCALE w_E scan      
QSTART
  Initial Q95 for RESCALE q_95 scan
QEND
  Final Q95 for RESCALE q_95 scan
