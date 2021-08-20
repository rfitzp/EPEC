# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage3 Lundquist numbers at rational surfaces

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = root['NEOCLASSICAL']['OUTPUTS']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
sk = ds['S']

ss = np.log10(np.asarray(sk))

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Lundquist Numbers")
psilow = psin[0] - 0.025

plt.xlim(psilow, 1.0)
plt.plot(psin, ss, 'r--', linewidth=0.4)
plt.plot(psin, ss, 'ko', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$log_{10}(S)$')

plt.show()
