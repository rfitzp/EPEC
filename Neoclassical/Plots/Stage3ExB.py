# -*-Python-*-
# Created by fitzpatrickr on 7 Aug 2021

# Script plots Stage3 ExB frequencies at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
web0 = ds['w_EB0']
web1 = ds['w_EB1']
web2 = ds['w_EB2']

fig = plt.figure(figsize=(12.0, 6.0))
fig.canvas.manager.set_window_title("FLUX: ExB Frequencies")
psilow = psin[0] - 0.025

plt.xlim(psilow, 1.0)
plt.plot(psin, web0, 'r--', linewidth=0.4)
plt.plot(psin, web0, 'ro', markersize=6)
plt.plot(psin, web1, 'g--', linewidth=0.4)
plt.plot(psin, web1, 'go', markersize=5)
plt.plot(psin, web2, 'b--', linewidth=0.4)
plt.plot(psin, web2, 'bo', markersize=5)
plt.axhline(y=0.0, color='y', linestyle='--')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\omega_{EB}(krad/s)$')

plt.show()
