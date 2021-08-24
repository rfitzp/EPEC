# -*-Python-*-
# Created by fitzpatrickr on 11 Aug 2021

# Script plots Stage5 island natural frequencies versus time

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = root['PHASE']['OUTPUTS']['Stage5']
ds = nc.Dataset(fn)
time = ds['time']
ome = ds['omega0']
mpol = ds['m_pol']

Time = np.asarray(time)
Ome = np.asarray(ome)
Mpol = np.asarray(mpol)

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("PHASE: Island Natural Frequencies (krad/s)")

ntime = Ome.shape[0]
nres = Ome.shape[1]

nrows = int(nres / 4)
if nrows * 4 < nres:
    nrows = nrows + 1

for n in range(nres):
    plt.subplot(nrows, 4, n + 1)
    plt.xlim(Time[0], Time[ntime - 1])
    ww = ome[:, n]
    plt.plot(Time, ww, color='blue', linewidth=2)
    plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
    plt.xlabel('$time(ms)$')
    strr = '$\\varpi_0 (m = ' + str(Mpol[n]) + ')$'
    plt.ylabel(strr)

plt.tight_layout(pad=0.5)

plt.show()
