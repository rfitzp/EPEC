import matplotlib.pyplot as plt
import numpy as np

D = 1.

D2 = D*D
D4 = D2*D2
D6 = D4*D2

def f1(q):
    return q*q/D2

def f2(q):
    return q**1.5

def f3(q):
    return q*D2

t1 = np.arange(0., D4, 1.e-3)
t2 = np.arange(D4, 2., 1.e-3)
t3 = np.arange(D4, 2., 1.e-3)

plt.hlines(y=D6, xmin=0, xmax=D4, linewidth=2, color='blue')
plt.plot(t1, f1(t1), linewidth=2, color='blue')
plt.plot(t2, f2(t2), linewidth=2, color='blue')
plt.plot(t3, f3(t3), linewidth=2, color='blue')

plt.text(0.4, 2.0, 'VR',  fontsize=20)
plt.text(1.5, 0.5, 'SC',  fontsize=20)
plt.text(0.2, 0.5, 'DR',  fontsize=20)
plt.text(1.75, 2.0, 'RI',  fontsize=20)

plt.text(1.29, 2.0, r'$P=Q_\ast^{3/2}$', fontsize=20, color='red')
plt.text(0.35, 1.03, '$P=D^6$', fontsize=20, color='red')
plt.text(0.741, 0.4, r'$P=Q_\ast^2 D^{-2}$', fontsize=20, color='red')
plt.text(1.525, 1.4, r'$P=Q_\ast D^2$', fontsize=20, color='red')

plt.xlabel(r'$Q_\ast/D^4$', fontsize=16)
plt.ylabel('$P/D^6$', fontsize=16)

plt.savefig("RegimeIII.pdf")

#plt.show()
