# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage2 safety factor profile

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

try:
    fn = root['OUTPUTS']['FLUX']['Stage2']
except:
    fn = '../Outputs/Stage2.nc'
ds = nc.Dataset(fn)

q = ds['q']
pn = ds['PsiN']
qr = ds['q_res']
pnr = ds['PsiN_res']
ppn = np.asarray(pn)
ppn = 1.0 - ppn

fig = plt.figure()
fig.canvas.manager.set_window_title("FLUX: Stage2 Safety Factor Profile")

plt.xlim((0.0, 1.0))
plt.plot(ppn, q, 'r--', linewidth=0.4)
plt.plot(pnr, qr, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$q$')

plt.show()
