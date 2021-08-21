# -*-Python-*-
# Created by fitzpatrickr on 4 Aug 2021

# Script plots Stage1 equilibirum

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['FLUX']['OUTPUTS']['Stage1']
ds = nc.Dataset(fn)

para = ds['Parameters']
rl = ds['RLIM']
zl = ds['ZLIM']
rlft = para[2]
rrgt = para[3]
zlow = para[4]
zhi = para[5]
rax = para[6]
zax = para[7]
pax = para[8]

aspect = (zhi - zlow) / (rrgt - rlft)

fn1 = root['FLUX']['OUTPUTS']['Stage2']
ds1 = nc.Dataset(fn1)
psir = ds1['PSI_R']
psiz = ds1['PSI_Z']
r = ds1['R']
z = ds1['Z']
rb = ds1['RBPTS']
zb = ds1['ZBPTS']

fig = plt.figure(figsize=(2.0 * 6.0 / aspect ** 0.5, 6.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title("FLUX: Stage2 Equilibrium Gradients")

XX, YY = np.meshgrid(r, z, indexing='ij')
ZR = np.asarray(psir)
zrmin = np.min(ZR)
zrmax = np.max(ZR)
levr = np.linspace(zrmin, zrmax, 40)
ZZ = np.asarray(psiz)
zzmin = np.min(ZZ)
zzmax = np.max(ZZ)
levz = np.linspace(zzmin, zzmax, 40)
zlevel = np.array([0.0])

plt.subplot(1, 2, 1)
plt.contour(XX, YY, ZR, levr)
plt.contour(XX, YY, ZR, zlevel, colors='red')
plt.plot(rax, zax, 'kx', markersize=4)
plt.plot(rb, zb, color='blue', linewidth=1)
plt.plot(rl, zl, color='black', linewidth=4)
plt.xlabel('$R/R_0$')
plt.ylabel("$Z/R_0$")
plt.title("$\\Psi_R$")

plt.subplot(1, 2, 2)
plt.contour(XX, YY, ZZ, levr)
plt.contour(XX, YY, ZZ, zlevel, colors='red')
plt.plot(rax, zax, 'kx', markersize=4)
plt.plot(rb, zb, color='blue', linewidth=1)
plt.plot(rl, zl, color='black', linewidth=4)
plt.xlabel('$R/R_0$')
plt.ylabel("$Z/R_0$")
plt.title("$\\Psi_Z$")

plt.tight_layout(pad=0.5)

plt.show()
