# -*-Python-*-
# Created by fitzpatrickr on 10 Aug 2021

# Script plots Stage4 vacuum island widths versus RMP coil phase

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = root['PHASE']['OUTPUTS']['Stage4']
ds = nc.Dataset(fn)
mpol = ds['m_pol']
pha = ds['phase']
wvac = ds['W_vacuum']

Mpol = np.asarray(mpol)
Phase = np.asarray(pha)
Wvac = np.asarray(wvac)

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("PHASE: Vacuum Island Widths")

nres = Wvac.shape[0]
nphi = Wvac.shape[1]

nrows = int(nres / 4)
if nrows * 4 < nres:
    nrows = nrows + 1

for n in range(nres):
    plt.subplot(nrows, 4, n + 1)
    plt.xlim(Phase[0], Phase[nphi - 1])
    ww = Wvac[n]
    plt.plot(Phase, ww, color='blue', linewidth=2)
    plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
    plt.xlabel('$\\phi_{rmp}/\\pi$')
    strr = '$W_{vac} (m = ' + str(Mpol[n]) + ')$'
    plt.ylabel(strr)

plt.tight_layout(pad=0.5)

plt.show()
