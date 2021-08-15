# -*-Python-*-
# Created by fitzpatrickr on 15 Aug 2021

# Script plots Stage1 raw plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

try:
    fn = root['OUTPUTS']['FLUX']['Stage1']
except:
    fn = '../Outputs/Stage1.nc'
ds = nc.Dataset(fn)

psin = ds['PSI_N']
p = ds['P']
pp = ds['Pp']
t = ds['T']
tp = ds['TTp']
q = ds['Q']
j = ds['J_phi']

fig = plt.figure(figsize=(8.0, 5.33))
fig.canvas.manager.set_window_title("FLUX: Stage1 Raw Profiles")

plt.subplot(2, 3, 1)
plt.xlim((0.0, 1.0))
plt.plot(psin, p)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$P$')

plt.subplot(2, 3, 2)
plt.xlim((0.0, 1.0))
plt.plot(psin, t)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$T$")

plt.subplot(2, 3, 3)
plt.xlim((0.0, 1.0))
plt.plot(psin, q)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$q$")

plt.subplot(2, 3, 4)
plt.xlim((0.0, 1.0))
plt.plot(psin, pp)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$P'$")

plt.subplot(2, 3, 5)
plt.xlim((0.0, 1.0))
plt.plot(psin, tp)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$TT'$")

plt.subplot(2, 3, 6)
plt.xlim((0.0, 1.0))
plt.plot(psin, j)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$J_\\phi$")

plt.tight_layout(pad=0.5)

plt.show()
