RESCALE NAMELIST
================

Short Description
-----------------

Description of RESCALE namelist

Keywords
--------

namelist

Namelist Variables
------------------

TYPE
   Rescaling type (1=density, 2=temperature, 3=size, 4=chi, 5=pressure, 6=wE, 7=q)
SCALE
   Rescaling factor (TYPE=1,2,3,4)
PSHIFT
   Pressure profile shift factor (kPa) (TYPE=5)
WSHIFT
   ExB frequency profile shift factor (krad/s) (TYPE=6
Q95
   Target safety-factor at 95% flux surface (TYPE=7)
OPOINT
   Flag for recalculating O-point position in rescaled equilibirum
XPOINT
   Flag for recalculating X-point position in rescaled equilibirum
