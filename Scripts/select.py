# Script to set up symbolic links for simulations of DIII-D discharge 145380

import sys
import string 
import os

narg = len (sys.argv) - 1
if narg == 0:
    time = input ("time (ms) ")
else:
    time = sys.argv[1]

basename   = "../ExperimentalData/DIIID/145380/gFiles/g145380.0"
null       = ""

filename   = null.join ([basename, time])

linkname   = "gFile"
space      = " "

remove     = "rm -f"
removelink = space.join ([remove, linkname])

link       = "ln -s"
createlink = space.join ([link, filename, linkname])

print     (removelink)
os.system (removelink)

print     (createlink)
os.system (createlink)

