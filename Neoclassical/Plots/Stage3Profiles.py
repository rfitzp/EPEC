# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage3 plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = root['NEOCLASSICAL']['OUTPUTS']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN']
n_e = ds['n_e']
T_e = ds['T_e']
n_i = ds['n_i']
T_i = ds['T_i']
Z_eff = ds['Z_eff']
n_n = ds['n_n']
w_E = ds['w_E']
w_t = ds['w_t']
w_p = ds['w_p']

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage3 Profiles")
psilow = 0.0

plt.subplot(3, 3, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, n_e)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$n_e(10^{19}\,m^{-3})$')

plt.subplot(3, 3, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, T_e)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$T_e(keV)$")

plt.subplot(3, 3, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, n_i)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$n_i(10^{19}\,m^{-3})$")

plt.subplot(3, 3, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, T_i)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$T_i(keV)$")

plt.subplot(3, 3, 5)
plt.xlim(psilow, 1.0)
plt.plot(psin, Z_eff)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$Z_{eff}$")

plt.subplot(3, 3, 6)
plt.xlim(psilow, 1.0)
plt.plot(psin, n_n)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$n_n((10^{19}\,m^{-3})$")

plt.subplot(3, 3, 7)
plt.xlim(psilow, 1.0)
plt.plot(psin, w_E)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_E(krad/s)$")

plt.subplot(3, 3, 8)
plt.xlim(psilow, 1.0)
plt.plot(psin, w_t)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\phi\,I}(krad/s)$")

plt.subplot(3, 3, 9)
plt.xlim(psilow, 1.0)
plt.plot(psin, w_p)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\omega_{\\theta\,I}(krad/s)$")

plt.tight_layout(pad=0.5)

plt.show()
