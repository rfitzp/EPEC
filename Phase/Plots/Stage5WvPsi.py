# -*-Python-*-
# Created by fitzpatrickr on 11 Aug 2021

# Script plots Stage5 vacuum island extends in PsiN and time

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn1 = root['FLUX']['OUTPUTS']['Stage2']
ds1 = nc.Dataset(fn1)
psin = ds1['PsiN_res']

fn = root['PHASE']['OUTPUTS']['Stage5']
ds = nc.Dataset(fn)
time = ds['time']
psim = ds['Psi_v-']
psip = ds['Psi_v+']

Time = np.asarray(time)
Psim = np.asarray(psim)
Psip = np.asarray(psip)

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("PHASE: Vacuum Island Extents")
psilow = psin[0] - 0.05

ntime = Psim.shape[0]
nres = Psim.shape[1]

for n in range(nres):
    plt.xlim(Time[0], Time[ntime - 1])
    plt.ylim(psilow, 1.0)
    pm = Psim[:, n]
    pp = Psip[:, n]
    px = pm.copy()
    for i in range(ntime):
        if i % 2 == 0:
            px[i] = pm[i]
        else:
            px[i] = pp[i]
    if n % 4 == 0:
        plt.plot(Time, pm, color='black', linewidth=0.5)
        plt.plot(Time, pp, color='black', linewidth=0.5)
        plt.plot(Time, px, color='black', linewidth=0.5)
    elif n % 4 == 1:
        plt.plot(Time, pm, color='red', linewidth=0.5)
        plt.plot(Time, pp, color='red', linewidth=0.5)
        plt.plot(Time, px, color='red', linewidth=0.5)
    elif n % 4 == 2:
        plt.plot(Time, pm, color='green', linewidth=0.5)
        plt.plot(Time, pp, color='green', linewidth=0.5)
        plt.plot(Time, px, color='green', linewidth=0.5)
    else:
        plt.plot(Time, pm, color='blue', linewidth=0.5)
        plt.plot(Time, pp, color='blue', linewidth=0.5)
        plt.plot(Time, px, color='blue', linewidth=0.5)

    plt.xlabel('$time(ms)$')
    strr = '$\\Psi_N$'
    plt.ylabel(strr)

plt.tight_layout(pad=0.5)

plt.show()
