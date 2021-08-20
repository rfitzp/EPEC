# -*-Python-*-
# Created by fitzpatrickr on 7 Aug 2021

# Script plots Stage3 lengths at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['NEOCLASSICAL']['OUTPUTS']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
rhos = ds['rho_s']
delt = ds['delta']
rhoe = ds['rho_theta_e']
rhoi = ds['rho_theta_i']

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Lengthscales")
psilow = psin[0] - 0.025

plt.subplot(2, 2, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, rhos, 'r--', linewidth=0.4)
plt.plot(psin, rhos, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\\rho_s(cm)$')

plt.subplot(2, 2, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, delt, 'r--', linewidth=0.4)
plt.plot(psin, delt, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\delta(cm)$")

plt.subplot(2, 2, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, rhoe, 'r--', linewidth=0.4)
plt.plot(psin, rhoe, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\rho_{\\theta\,e}(cm)$")

plt.subplot(2, 2, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, rhoi, 'r--', linewidth=0.4)
plt.plot(psin, rhoi, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\rho_{\\theta\,i}(cm)$")

plt.tight_layout(pad=0.5)

plt.show()
