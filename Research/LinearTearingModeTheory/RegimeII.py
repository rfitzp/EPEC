import matplotlib.pyplot as plt
import numpy as np

D = 1.2

D2 = D*D
DM2 = 1./D2
D3 = D2*D
D4 = D2*D2
D6 = D4*D2

def f1(q):
    return q*q/D2

def f2(q):
    return D2/q/q

def f3(q):
    return q**3.

def f4(q):
    return 1./q**3.

def f5(q):
    return D4/q

t1 = np.arange(0., D, 1.e-3)
t2 = np.arange(DM2, D, 1.e-3)
t3 = np.arange(D, 4.**(1./3.), 1.e-3)
t4 = np.arange(4.**(-1./3.), DM2, 1.e-3)
t5 = np.arange(DM2, D, 1.e-3)

plt.hlines(y=D6, xmin=0, xmax=DM2, linewidth=2, color='blue')
plt.vlines(x=D, ymin=0, ymax=D3, linewidth=2, color='blue')
plt.plot(t1, f1(t1), linewidth=2, color='blue')
plt.plot(t2, f2(t2), linewidth=2, color='blue')
plt.plot(t3, f3(t3), linewidth=2, color='blue')
plt.plot(t4, f4(t4), linewidth=2, color='blue')
plt.plot(t5, f5(t5), linewidth=2, color='blue')

plt.text(0.3, 3.5, 'VR',  fontsize=20)
plt.text(1.1, 3.5, 'VI',  fontsize=20)
plt.text(1.4, 1.0, 'I',  fontsize=20)
plt.text(1.025, 1.5, 'DI',  fontsize=20)
plt.text(0.4, 1.0, 'DR',  fontsize=20)
plt.text(0.9, 0.2, 'SC',  fontsize=20)

plt.text(1.34, 2.25, '$P=Q^3$', fontsize=20, color='red')
plt.text(0.665, 3.5, '$P=Q^{-3}$', fontsize=20, color='red')
plt.text(0.58, 1.5, '$P=D^2\!Q^{-2}$', fontsize=20, color='red')
plt.text(1.21, 0.5, '$Q=D$', fontsize=20, color='red')
plt.text(0.25, 3.04, '$P=D^6$', fontsize=20, color='red')
plt.text(0.38, 0.3, '$P=Q^2 \!D^{-2}$', fontsize=20, color='red')
plt.text(0.91, 2.3, '$P=D^4\!Q^{-1}$', fontsize=20, color='red')

plt.xlabel('$Q$', fontsize=16)
plt.ylabel('$P$', fontsize=16)

plt.savefig("RegimeII.pdf")

#plt.show()
