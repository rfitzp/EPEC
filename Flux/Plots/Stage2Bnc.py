# -*-Python-*-
# Created by fitzpatrickr on 9 Aug 2021

# Script plots Stage2 B(Theta) at rational surfaces

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import math

fn = root['FLUX']['OUTPUTS']['Stage2']
ds = nc.Dataset(fn)
bnc = ds['B_nc']
cnc = ds['C_nc']

Bnc = np.asarray(bnc)
Cnc = np.asarray(cnc)

fig = plt.figure(figsize=(12.0, 6.0))
fig.canvas.manager.set_window_title("FLUX: Stage2 B(Theta)")

nres = Bnc.shape[0]
nTheta = Bnc.shape[1]
Theta = np.linspace(0.0, 2.0, nTheta)

plt.subplot(1, 2, 1)
plt.xlim(0.0, 2.0)
for i in range(nres):
    ttt = Bnc[i]
    plt.plot(Theta, ttt, color='blue', linewidth=0.5)
# plt.axhline (1., color='red', linewidth=0.5, linestyle='dotted')
plt.axvline(1.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Theta/\\pi$')
plt.ylabel("$B$")

plt.subplot(1, 2, 2)
plt.xlim(0.0, 2.0)
for i in range(nres):
    ttt = Cnc[i]
    plt.plot(Theta, ttt, color='blue', linewidth=0.5)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.axvline(1.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Theta/\\pi$')
plt.ylabel("$\\partial B/\\partial\\Theta$")

plt.show()
