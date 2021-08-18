# -*-Python-*-
# Created by fitzpatrickr on 13 Aug 2021

# Script plots Stage3 transport factors at rational surfaces

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

try:
    fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
except:
    fn = '../Outputs/Stage3.nc'
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
fac1 = ds['Factor_1']
fac2 = ds['Factor_2']
fac3 = ds['Factor_3']
fac4 = ds['Factor_4']
fac5 = ds['Factor_5']
fac6 = ds['Factor_6']
fac7 = ds['Factor_7']
fac8 = ds['Factor_8']
fac9 = ds['Factor_9']
fac10 = ds['Factor_10']
fac11 = ds['Factor_11']
fac12 = ds['Factor_12']

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Transport Factors (10^19 m^-3 keV)")
psilow = psin[0] - 0.025

plt.subplot(3, 4, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac1, 'r--', linewidth=0.4)
plt.plot(psin, fac1, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 1')

plt.subplot(3, 4, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac2, 'r--', linewidth=0.4)
plt.plot(psin, fac2, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 2')

plt.subplot(3, 4, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac3, 'r--', linewidth=0.4)
plt.plot(psin, fac3, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 3')

plt.subplot(3, 4, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac4, 'r--', linewidth=0.4)
plt.plot(psin, fac4, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 4')

plt.subplot(3, 4, 5)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac5, 'r--', linewidth=0.4)
plt.plot(psin, fac5, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 5')

plt.subplot(3, 4, 6)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac6, 'r--', linewidth=0.4)
plt.plot(psin, fac6, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 6')

plt.subplot(3, 4, 7)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac7, 'r--', linewidth=0.4)
plt.plot(psin, fac7, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 7')

plt.subplot(3, 4, 8)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac8, 'r--', linewidth=0.4)
plt.plot(psin, fac8, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 8')

plt.subplot(3, 4, 9)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac9, 'r--', linewidth=0.4)
plt.plot(psin, fac9, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 9')

plt.subplot(3, 4, 10)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac10, 'r--', linewidth=0.4)
plt.plot(psin, fac10, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 10')

plt.subplot(3, 4, 11)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac11, 'r--', linewidth=0.4)
plt.plot(psin, fac11, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 11')

plt.subplot(3, 4, 12)
plt.xlim(psilow, 1.0)
plt.plot(psin, fac12, 'r--', linewidth=0.4)
plt.plot(psin, fac12, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('Factor 12')

plt.tight_layout(pad=0.5)

plt.show()
