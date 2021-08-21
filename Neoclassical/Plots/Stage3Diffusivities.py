# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage3 diffusivity profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = root['NEOCLASSICAL']['OUTPUTS']['Stage3']
ds = nc.Dataset(fn)

psin = ds['PsiN']
chie = ds['chi_e']
chii = ds['chi_i']
chip = ds['chi_p']
chin = ds['chi_n']

fig = plt.figure(figsize=(8.0, 5.33))
fig.canvas.manager.set_window_title("FLUX: Stage3 Diffusivities")
psilow = 0.0

plt.subplot(2, 2, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, chie)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\chi_e(m^2/s)$')

plt.subplot(2, 2, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, chii)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\chi_i(m^2/s)$")

plt.subplot(2, 2, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, chip)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\chi_\phi(m^2/s)$")

plt.subplot(2, 2, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, chin)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$D_\perp(m^2/s)$")

plt.tight_layout(pad=0.5)

plt.show()
