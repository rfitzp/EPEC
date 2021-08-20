# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage2 D_R profile

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = root['FLUX']['OUTPUTS']['Stage2']
ds = nc.Dataset(fn)

dr = ds['DR_res']
pnr = ds['PsiN_res']

fig = plt.figure()
fig.canvas.manager.set_window_title("FLUX: Stage2 D_R Profile")
psilow = pnr[0] - 0.025

plt.xlim(psilow, 1.0)
plt.plot(pnr, dr, 'r--', linewidth=0.4)
plt.plot(pnr, dr, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$D_R$')

plt.show()
