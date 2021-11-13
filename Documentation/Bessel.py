import scipy.special as sc
import matplotlib.pyplot as plt
import numpy as np
import math

N = 1000
eps = 0.01

j0p = sc.jn_zeros(0, N)
j1p = sc.jn_zeros(1, N)

def Sum1(x):

    sum = 0.
    for i in range(0, N):
        sum = sum + sc.jn(0, j0p[i]*x) * sc.jn(0, j0p[i]*x) /sc.jn(1, j0p[i])/sc.jn(1, j0p[i]) /j0p[i]/j0p[i]

    return sum

sum1 = np.vectorize(Sum1)

def Log(x):

    return 0.5*math.log(1./x)

logg = np.vectorize(Log)

def Sum2(x):

    sum = 0.
    for i in range(0, N):
        sum = sum + sc.jn(1, j1p[i]*x) * sc.jn(1, j1p[i]*x) /sc.jn(2, j1p[i])/sc.jn(2, j1p[i]) /(1. + eps*j1p[i]*j1p[i])

    return sum*eps**0.5

sum2 = np.vectorize(Sum2)

def Exp(x):

    return 0.25/x

expp = np.vectorize(Exp)

x1 = np.arange (0.01, 1., 1.e-3)
np.append(x1, 1.)

plt.plot(x1, sum1(x1), linewidth=2, color='blue')
plt.plot(x1, logg(x1), linewidth=2, linestyle='dashed', color='red')

x2 = np.arange (0.01, 1., 1.e-3)
np.append(x2, 1.)

plt.xlim(0.,1.)
plt.ylim(0.,3.)
plt.xlabel(r'$r_s/a$', fontsize=16)
plt.ylabel(r'$F_\theta,\, F_\varphi$', fontsize=16)

eps = 1.e-3
#plt.plot(x2, sum2(x2), linewidth=2, color='green')

eps = 1.e-4
plt.plot(x2, sum2(x2), linewidth=2, color='green')

eps = 1.e-5
plt.plot(x2, expp(x2), linewidth=2, linestyle='dashed', color='red')

#plt.savefig("Bessel.pdf")

plt.show()

