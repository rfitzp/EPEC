# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage1 plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = root['FLUX']['OUTPUTS']['Stage1']
ds = nc.Dataset(fn)

psin = ds['PSI_N']
p = ds['p']
pp = ds['pp']
t = ds['t']
tp = ds['ttp']
q = ds['q']
j = ds['j_phi']

fig = plt.figure(figsize=(8.0, 5.33))
fig.canvas.manager.set_window_title("FLUX: Stage1 Profiles")

plt.subplot(2, 3, 1)
plt.xlim((0.0, 1.0))
plt.plot(psin, p)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$P$')

plt.subplot(2, 3, 2)
plt.xlim((0.0, 1.0))
plt.plot(psin, t)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$g$")

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
plt.ylabel("$gg'$")

plt.subplot(2, 3, 6)
plt.xlim((0.0, 1.0))
plt.plot(psin, j)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$J_\\phi$")

plt.tight_layout(pad=0.5)

plt.show()
