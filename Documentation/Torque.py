import matplotlib.pyplot as plt
import numpy as np

zeta = 5.e-2
I    = 1.

def Tvs(v):

    return 1. - v

tvs = np.vectorize(Tvs)

def Tem(v):

    return 0.25 * I*I * v /(zeta*zeta + v*v)

tem = np.vectorize(Tem)

v1 = np.arange(0., 1., 1.e-3)
np.append(v1, 1.)

plt.xlim(0.,1.)
plt.ylim(0.,2.6)
plt.xlabel(r'$v$', fontsize=16)
plt.ylabel(r'$T_{\rm VS},\, T_{\rm EM}$', fontsize=16)

I = 0.8
plt.plot(v1, tem(v1), linewidth=2, color='blue')

plt.plot(v1, tvs(v1), linewidth=2, color='red')

I = 1.01
plt.plot(v1, tem(v1), linewidth=2, linestyle='dashed',color='blue')

plt.text(0.72, 1., 'shielded', fontsize=20, color='black')
plt.arrow(0.8, 0.97, 0., -0.72)

plt.text(0.175, 1.5, 'unshielded', fontsize=20, color='black')
plt.arrow(0.25, 1.47, -0.22, -0.46)

plt.text(0.41, 0.8, 'bifurcation', fontsize=20, color='black')
plt.arrow(0.5, 0.77, 0., -0.22)

plt.savefig("Torque.pdf")

plt.show()

