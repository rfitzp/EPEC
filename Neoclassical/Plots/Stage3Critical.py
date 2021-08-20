# -*-Python-*-
# Created by fitzpatrickr on 7 Aug 2021

# Script plots Stage3 critical island widths at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
wcre = ds['W_crit_Te']
wcri = ds['W_crit_Ti']
wcrn = ds['W_crit_ne']

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Critical Island Widths")
psilow = psin[0] - 0.025

plt.subplot(3, 1, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, wcre, 'r--', linewidth=0.4)
plt.plot(psin, wcre, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$W_{crit\,T_e}(cm)$')

plt.subplot(3, 1, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, wcri, 'r--', linewidth=0.4)
plt.plot(psin, wcri, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$W_{crit\,T_i}(cm)$")

plt.subplot(3, 1, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, wcrn, 'r--', linewidth=0.4)
plt.plot(psin, wcrn, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$W_{crit\,n_e}(cm)$")

plt.tight_layout(pad=0.5)

plt.show()
