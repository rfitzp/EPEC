# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage3 natural frequencies at rational surfaces

import netCDF4 as nc
import matplotlib.pyplot as plt

try:
    fn = root['OUTPUTS']['NEOCLASSICAL']['Stage3']
except:
    fn = '../Outputs/Stage3.nc'
ds = nc.Dataset(fn)

psin = ds['PsiN_k']
wl = ds['w_linear']
wn = ds['w_nonlinear']
we = ds['w_EB']

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Natural Frequencies")
psilow = psin[0] - 0.025

plt.xlim(psilow, 1.0)
plt.plot(psin, wl, 'r--', linewidth=0.4)
plt.plot(psin, wl, 'ro', markersize=5)
plt.plot(psin, wn, 'b--', linewidth=0.4)
plt.plot(psin, wn, 'bo', markersize=5)
plt.plot(psin, we, 'k--', linewidth=0.4)
plt.plot(psin, we, 'ko', markersize=5)
plt.axhline(y=0.0, color='y', linestyle='--')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\\varpi(krad/s)$')

plt.show()
