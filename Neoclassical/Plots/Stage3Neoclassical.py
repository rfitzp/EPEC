# -*-Python-*-
# Created by fitzpatrickr on 14 Aug 2021

# Script plots Stage3 neoclassical frequencies at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
wnci = ds['w_nc_i']
wncI0 = ds['w_nc_I0']
wncI1 = ds['w_nc_I1']
wncI2 = ds['w_nc_I2']
wnce = ds['w_nc_e']
wthi = ds['w_th_i']
wthI = ds['w_theta_I']
wthI0 = ds['w_pnc_I0']
wthI0 = ds['w_pnc_I0']
wthI1 = ds['w_pnc_I1']
wthI2 = ds['w_pnc_I2']
wpar = ds['w_para']
wast = ds['w_ast_i']
w = -np.asarray(wast)

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Neoclassical Frequencies (krad/s)")
psilow = psin[0] - 0.025

plt.subplot(2, 3, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, wnci, 'r--', linewidth=0.4)
plt.plot(psin, wnci, 'bo', markersize=5)
plt.plot(psin, w, 'r--', linewidth=0.4)
plt.plot(psin, w, 'ro', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\\omega_{nc\,i}$')

plt.subplot(2, 3, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, wncI0, 'r--', linewidth=0.4)
plt.plot(psin, wncI0, 'ro', markersize=5)
plt.plot(psin, wncI1, 'r--', linewidth=0.4)
plt.plot(psin, wncI1, 'go', markersize=5)
plt.plot(psin, wncI2, 'r--', linewidth=0.4)
plt.plot(psin, wncI2, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{nc\,I}$")

plt.subplot(2, 3, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, wnce, 'r--', linewidth=0.4)
plt.plot(psin, wnce, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{nc\,e}$")

plt.subplot(2, 3, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, wthi, 'r--', linewidth=0.4)
plt.plot(psin, wthi, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\\theta\,i}$")

plt.subplot(2, 3, 5)
plt.xlim(psilow, 1.0)
plt.plot(psin, wthI, 'r--', linewidth=0.4)
plt.plot(psin, wthI, 'ko', markersize=5)
plt.plot(psin, wthI0, 'r--', linewidth=0.4)
plt.plot(psin, wthI0, 'ro', markersize=5)
plt.plot(psin, wthI1, 'r--', linewidth=0.4)
plt.plot(psin, wthI1, 'go', markersize=5)
plt.plot(psin, wthI2, 'r--', linewidth=0.4)
plt.plot(psin, wthI2, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\\theta\,I}$")

plt.subplot(2, 3, 6)
plt.xlim(psilow, 1.0)
plt.plot(psin, wpar, 'r--', linewidth=0.4)
plt.plot(psin, wpar, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\\parallel}$")

plt.tight_layout(pad=0.5)

plt.show()
