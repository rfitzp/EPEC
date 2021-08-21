# -*-Python-*-
# Created by fitzpatrickr on 18 Aug 2021

# Script plots Stage3 right-hand-sides of Rutherford equations

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = root['NEOCLASSICAL']['OUTPUTS']['Stage3']
ds = nc.Dataset(fn)
mpol = ds['m_k']
w = ds['w_psi']
rhs = ds['rhs']

Mpol = np.asarray(mpol)
W = np.asarray(w)
RHS = np.asarray(rhs)

fig = plt.figure(figsize=(12.0, 8.0))
fig.canvas.manager.set_window_title("NEOCLASSICAL: Right-Hand Sides of Rutherford Equations")

nres = RHS.shape[0]
nw = RHS.shape[1]

nrows = int(nres / 4)
if nrows * 4 < nres:
    nrows = nrows + 1

for n in range(nres):
    plt.subplot(nrows, 4, n + 1)
    plt.xlim(0.0, W[nw - 1])
    ww = rhs[n]
    plt.plot(W, ww, color='blue', linewidth=2)
    plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
    plt.xlabel('$W$')
    strr = 'RHS $(m = ' + str(Mpol[n]) + ')$'
    plt.ylabel(strr)

plt.tight_layout(pad=0.5)

plt.show()
