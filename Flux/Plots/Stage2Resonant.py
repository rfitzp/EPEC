# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage2 resonant quantities

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = root['FLUX']['OUTPUTS']['Stage2']
ds = nc.Dataset(fn)

psin = ds['PsiN_res']
q = ds['q_res']
s = ds['s_res']
g = ds['g_res']
gm = ds['gamma_res']
fc = ds['fc_res']
aj = ds['ajj_res']
Q = ds['Q_res']
A1 = ds['A1_res']
A2 = ds['A2_res']

fig = plt.figure(figsize=(8.0, 8.0))
fig.canvas.manager.set_window_title("FLUX: Stage2 Resonant Quantities")
psilow = psin[0] - 0.025

plt.subplot(3, 3, 1)
plt.xlim(psilow, 1.0)
plt.plot(psin, q, 'r--', linewidth=0.4)
plt.plot(psin, q, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$q$')

plt.subplot(3, 3, 2)
plt.xlim(psilow, 1.0)
plt.plot(psin, s, 'r--', linewidth=0.4)
plt.plot(psin, s, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$s$")

plt.subplot(3, 3, 3)
plt.xlim(psilow, 1.0)
plt.plot(psin, g, 'r--', linewidth=0.4)
plt.plot(psin, g, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$g$")

plt.subplot(3, 3, 4)
plt.xlim(psilow, 1.0)
plt.plot(psin, gm, 'r--', linewidth=0.4)
plt.plot(psin, gm, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\gamma$")

plt.subplot(3, 3, 5)
plt.xlim(psilow, 1.0)
plt.plot(psin, fc, 'r--', linewidth=0.4)
plt.plot(psin, fc, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$f_c$")

plt.subplot(3, 3, 6)
plt.xlim(psilow, 1.0)
plt.plot(psin, aj, 'r--', linewidth=0.4)
plt.plot(psin, aj, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$a_{jj}$")

plt.subplot(3, 3, 7)
plt.xlim(psilow, 1.0)
plt.plot(psin, Q, 'r--', linewidth=0.4)
plt.plot(psin, Q, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$Q$")

plt.subplot(3, 3, 8)
plt.xlim(psilow, 1.0)
plt.plot(psin, A1, 'r--', linewidth=0.4)
plt.plot(psin, A1, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$A_1$")

plt.subplot(3, 3, 9)
plt.xlim(psilow, 1.0)
plt.plot(psin, A2, 'r--', linewidth=0.4)
plt.plot(psin, A2, 'bo', markersize=3)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$A_2$")

plt.tight_layout(pad=0.5)

plt.show()
