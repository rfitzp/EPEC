# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage1 plasma equilibrium

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = root['FLUX']['OUTPUTS']['Stage1']
ds = nc.Dataset(fn)

para = ds['Parameters']
psi = ds['PSI']
r = ds['R']
z = ds['Z']
rb = ds['RBOUND']
zb = ds['ZBOUND']
rl = ds['RLIM']
zl = ds['ZLIM']
rlft = para[2]
rrgt = para[3]
zlow = para[4]
zhi = para[5]
rax = para[6]
zax = para[7]
pax = para[8]
pab = para[9]

aspect = (zhi - zlow) / (rrgt - rlft)

fig = plt.figure(figsize=(6.0 / aspect ** 0.5, 6.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title("FLUX: Stage1 Equilibrium")

XX, YY = np.meshgrid(r, z, indexing='ij')
ZZ = np.asarray(psi)
levels = np.linspace(pax, pab, 40)

plt.contour(XX, YY, ZZ, levels)

plt.plot(rax, zax, 'ro', markersize=2)
plt.plot(rb, zb, color='blue', linewidth=1)
plt.plot(rl, zl, color='black', linewidth=4)
plt.xlabel('$R/R_0$')
plt.ylabel("$Z/R_0$")

plt.show()
