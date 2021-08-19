# -*-Python-*-
# Created by fitzpatrickr on 11 Aug 2021

# Script plots Stage5 vacuum island widths versus time

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

try:
    fn = root['OUTPUTS']['PHASE']['Stage4']
except:
    fn = '../Outputs/Stage4.nc'
ds = nc.Dataset(fn)
mpol = ds['m_pol']

Mpol = np.asarray(mpol)

try:
    fn1 = root['OUTPUTS']['PHASE']['Stage5']
except:
    fn1 = '../Outputs/Stage5.nc'
ds1 = nc.Dataset(fn1)
time = ds1['time']
wv = ds1['W_vac']
w = ds1['W']

Time = np.asarray(time)
Wv = np.asarray(wv)
W = np.asarray(w)

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("PHASE: Vacuum Island Widths")

ntime = Wv.shape[0]
nres = Wv.shape[1]

nrows = int(nres / 4)
if nrows * 4 < nres:
    nrows = nrows + 1

for n in range(nres):
    plt.subplot(nrows, 4, n + 1)
    plt.xlim(Time[0], Time[ntime - 1])
    ww = Wv[:, n]
    w1 = W[:, n]
    plt.plot(Time, ww, color='blue', linewidth=2)
    plt.plot(Time, w1, color='red', linewidth=2)
    plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
    plt.xlabel('$time(ms)$')
    strr = '$W_{vac} (m = ' + str(Mpol[n]) + ')$'
    plt.ylabel(strr)

plt.tight_layout(pad=0.5)

plt.show()
