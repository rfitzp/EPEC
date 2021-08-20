# -*-Python-*-
# Created by fitzpatrickr on 11 Aug 2021

# Script plots Stage5 RMP waveform

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = root['OUTPUTS']['PHASE']['Stage5']
ds = nc.Dataset(fn)

time = ds['time']
irmp = ds['Irmp']
prmp = ds['Prmp']

Time = np.asarray(time)

ntim = Time.shape[0]

fig = plt.figure(figsize=(8.0, 4.0))
fig.canvas.manager.set_window_title("PHASE: Stage5 RMP Waveform")

plt.subplot(1, 2, 1)
plt.xlim((Time[0], Time[ntim - 1]))
plt.plot(time, irmp, color='blue', linewidth=1)
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.xlabel('$time(ms)$')
plt.ylabel('$I_{rmp}(kA)$')

plt.subplot(1, 2, 2)
plt.xlim((Time[0], Time[ntim - 1]))
plt.ylim((-1.0, 1.0))
plt.axhline(0.0, color='red', linewidth=0.5, linestyle='dotted')
plt.plot(time, prmp, color='blue', linewidth=1)
plt.xlabel('$time(ms)$')
plt.ylabel("$\\phi_{rmp}/\\pi$")

plt.tight_layout(pad=0.5)

plt.show()
