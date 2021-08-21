# -*-Python-*-
# Created by fitzpatrickr on 9 Aug 2021

# Script plots Stage2 theta(Theta) at rational surfaces

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import math

fn = root['FLUX']['OUTPUTS']['Stage2']
ds = nc.Dataset(fn)
th = ds['theta']

theta = np.asarray(th) / math.pi

fig = plt.figure(figsize=(6.0, 6.0))
fig.canvas.manager.set_window_title("FLUX: Stage2 theta(Theta)")

nres = theta.shape[0]
nTheta = theta.shape[1]
Theta = np.linspace(0.0, 2.0, nTheta)

plt.xlim(0.0, 2.0)
plt.ylim(0.0, 2.0)
for i in range(nres):
    ttt = theta[i]
    plt.plot(Theta, ttt, color='blue', linewidth=0.5)
plt.axhline(1.0, color='red', linewidth=0.5, linestyle='dotted')
plt.axvline(1.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Theta/\\pi$')
plt.ylabel("$\\theta/\\pi$")

plt.show()
