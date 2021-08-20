# -*-Python-*-
# Created by fitzpatrickr on 8 Aug 2021

# Script plots Stage1 straight coordinate sdstem

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

fn = root['OUTPUTS']['FLUX']['Stage1']
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

fn1 = root['OUTPUTS']['FLUX']['Stage2']
ds1 = nc.Dataset(fn1)
rb = ds1['RBPTS']
zb = ds1['ZBPTS']
rst = ds1['R_st']
zst = ds1['Z_st']

R = np.asarray(rst)
Z = np.asarray(zst)

fig = plt.figure(figsize=(6.0 / aspect ** 0.5, 6.0 * aspect ** 0.5))
fig.canvas.manager.set_window_title("FLUX: Stage2 Straight Coordinates")

plt.subplot(1, 1, 1)
nth = R.shape[1]
for j in range(0, nth, 1):
    rr1 = R[:, j]
    zz1 = Z[:, j]
    plt.plot(rr1, zz1, color='black', linewidth=0.5)
nres = R.shape[0]
for i in range(nres):
    rr = R[i]
    zz = Z[i]
    plt.plot(rr, zz, color='red', linewidth=1.0)
plt.plot(rax, zax, 'kx', markersize=4)
plt.plot(rb, zb, color='blue', linewidth=1)
plt.xlabel('$R/R_0$')
plt.ylabel("$Z/R_0$")

plt.show()