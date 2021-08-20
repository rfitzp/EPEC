# -*-Python-*-
# Created by fitzpatrickr on 11 Aug 2021

# Script plots Stage5 island phases versus time

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = root['OUTPUTS']['PHASE']['Stage4']
ds = nc.Dataset(fn)
mpol = ds['m_pol']

Mpol = np.asarray(mpol)

fn1 = root['OUTPUTS']['PHASE']['Stage5']
ds1 = nc.Dataset(fn1)
time = ds1['time']
phi = ds1['phi']

Time = np.asarray(time)
Phi = np.asarray(phi)

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("PHASE: Island Phases")

ntime = Phi.shape[0]
nres = Phi.shape[1]

nrows = int(nres / 4)
if nrows * 4 < nres:
    nrows = nrows + 1

for n in range(nres):
    plt.subplot(nrows, 4, n + 1)
    plt.xlim(Time[0], Time[ntime - 1])
    plt.ylim(-1.0, 1.0)
    ww = phi[:, n]
    plt.plot(Time, ww, color='blue', linewidth=2)
    plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
    plt.xlabel('$time(ms)$')
    strr = '$\phi/\\pi (m = ' + str(Mpol[n]) + ')$'
    plt.ylabel(strr)

plt.tight_layout(pad=0.5)

plt.show()
