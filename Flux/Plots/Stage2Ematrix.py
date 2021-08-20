# -*-Python-*-
# Created by fitzpatrickr on 9 Aug 2021

# Script plots Stage2 E-matrix

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import math

fn = root['OUTPUTS']['FLUX']['Stage2']
ds = nc.Dataset(fn)
er = ds['E_real']
ei = ds['E_imag']
mp = np.asarray(ds['m_pol'])

fig = plt.figure(figsize=(12.0, 6.0))
fig.canvas.manager.set_window_title("FLUX: Stage2 E-matrix")

Ereal = np.asarray(er)
Eimag = np.asarray(ei)

nres = mp.shape[0]
for i in range(nres):
    for j in range(nres):
        Ereal[i][j] = Ereal[i][j] / 2.0 / math.sqrt(float(mp[i] * mp[j]))
        Eimag[i][j] = Eimag[i][j] / 2.0 / math.sqrt(float(mp[i] * mp[j]))

plt.subplot(1, 2, 1)
plt.matshow(Ereal, fignum=0)
plt.colorbar()
plt.title('$E_{real}\, /2\,(mm\')^{1/2}$', pad=25)

plt.subplot(1, 2, 2)
plt.matshow(Eimag, fignum=0)
plt.colorbar()
plt.title('$E_{imag}\, /2\,(mm\')^{1/2}$', pad=25)

plt.show()
