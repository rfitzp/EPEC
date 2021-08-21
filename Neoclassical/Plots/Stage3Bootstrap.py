# -*-Python-*-
# Created by fitzpatrickr on 9 Aug 2021

# Script plots Stage3 bootstrap parameters at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['NEOCLASSICAL']['OUTPUTS']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
alpe = ds['alpha_b_e']
alpi = ds['alpha_b_i']
alpc = ds['alpha_c']
alpp = ds['alpha_p']
alpbc = np.asarray(alpe) + np.asarray(alpc)

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Lengthscales")
psilow = psin[0] - 0.025

plt.subplot(2, 2, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, alpe, 'r--', linewidth=0.4)
plt.plot(psin, alpe, 'bo', markersize=5)
plt.plot(psin, alpbc, 'b--', linewidth=0.4)
plt.plot(psin, alpbc, 'ro', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\\alpha_{b\,e}$')

plt.subplot(2, 2, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, alpi, 'r--', linewidth=0.4)
plt.plot(psin, alpi, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\alpha_{b\,i}$")

plt.subplot(2, 2, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, alpc, 'r--', linewidth=0.4)
plt.plot(psin, alpc, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\alpha_c$")

plt.subplot(2, 2, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, alpp, 'r--', linewidth=0.4)
plt.plot(psin, alpp, 'bo', markersize=5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\alpha_p$")

plt.tight_layout(pad=0.5)

plt.show()
