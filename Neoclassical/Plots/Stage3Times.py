# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage3 timescales at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
tauh = ds['tau_H']
taur = ds['tau_R']
taum = ds['tau_M']
taut = ds['tau_th']
taux = ds['tau_cx']

tauH = np.log10(np.asarray(tauh))
tauR = np.log10(np.asarray(taur))
tauM = np.log10(np.asarray(taum))
tauT = np.log10(np.asarray(taut))
tauX = np.log10(np.asarray(taux))

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Timescales")
psilow = psin[0] - 0.025

plt.subplot(3, 2, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, tauH, 'r--', linewidth=0.4)
plt.plot(psin, tauH, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$log_{10}(\\tau_H (s))$')

plt.subplot(3, 2, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, tauR, 'r--', linewidth=0.4)
plt.plot(psin, tauR, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$log_{10}(\\tau_R (s))$")

plt.subplot(3, 2, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, tauM, 'r--', linewidth=0.4)
plt.plot(psin, tauM, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$log_{10}(\\tau_M (s))$")

plt.subplot(3, 2, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, tauT, 'r--', linewidth=0.4)
plt.plot(psin, tauT, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$log_{10}(\\tau_\\theta (s))$")

plt.subplot(3, 2, 5)
plt.xlim(psilow, 1.0)
plt.plot(psin, tauX, 'r--', linewidth=0.4)
plt.plot(psin, tauX, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\log_{10}(\\tau_{cx} (s))$")

plt.tight_layout(pad=0.5)

plt.show()
