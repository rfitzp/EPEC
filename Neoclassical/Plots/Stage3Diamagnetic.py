# -*-Python-*-
# Created by fitzpatrickr on 7 Aug 2021

# Script plots Stage3 diamagnetic frequencies at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
wase = ds['w_ast_e']
wasi = ds['w_ast_i']
wasI = ds['w_ast_I']

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Diamagnetic Frequencies")
psilow = psin[0] - 0.025

plt.subplot(3, 1, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, wase, 'r--', linewidth=0.4)
plt.plot(psin, wase, 'bo', markersize=5)
plt.axhline(y=0.0, color='y', linestyle='--')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\omega_{\\ast\,e}(krad/s)$')

plt.subplot(3, 1, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, wasi, 'r--', linewidth=0.4)
plt.plot(psin, wasi, 'bo', markersize=5)
plt.axhline(y=0.0, color='y', linestyle='--')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\\ast\,i}(krad/s)$")

plt.subplot(3, 1, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, wasI, 'r--', linewidth=0.4)
plt.plot(psin, wasI, 'bo', markersize=5)
plt.axhline(y=0.0, color='y', linestyle='--')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\\ast\,I}(krad/s)$")

plt.tight_layout(pad=0.5)

plt.show()
