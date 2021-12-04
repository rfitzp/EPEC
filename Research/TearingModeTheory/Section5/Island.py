import matplotlib.pyplot as plt
import numpy as np
import math

F = 0.

def xp(t):
    return math.sqrt(F - math.cos(t*math.pi))/2.**0.5/2.

Xplus = np.vectorize(xp)

def xm(t):
    return - math.sqrt(F - math.cos(t*math.pi))/2.**0.5/2.

Xminus = np.vectorize(xm)

def xxp(t):
    return  math.sqrt(-math.cos(t*math.pi)) /2./2.**0.5

XXplus = np.vectorize(xxp)

def xxm(t):
    return  - math.sqrt(-math.cos(t*math.pi)) /2./2.**0.5

XXminus= np.vectorize(xxm)

t1 = np.arange(0., 6., 1.e-3)
t2 = np.arange(0.5000000001, 1.49999999999, 1.e-4)
t3 = t2 + 2.
t4 = t3 + 2.

F = 1.
plt.plot(t1, Xplus(t1), linewidth=2, color='red')
plt.plot(t1, Xminus(t1), linewidth=2, color='red')

F = 2.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 3.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 4.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 5.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 6.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 7.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 8.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 9.0
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 10.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 11.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 12.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 13.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 14.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 15.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 16.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 17.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 18.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 19.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')

F = 20.
plt.plot(t1, Xplus(t1), linewidth=2, color='blue')
plt.plot(t1, Xminus(t1), linewidth=2, color='blue')


plt.plot(t2, XXplus(t2), linewidth=2, color='blue')
plt.plot(t2, XXminus(t2), linewidth=2, color='blue')
plt.plot(t3, XXplus(t2), linewidth=2, color='blue')
plt.plot(t3, XXminus(t2), linewidth=2, color='blue')
plt.plot(t4, XXplus(t2), linewidth=2, color='blue')
plt.plot(t4, XXminus(t2), linewidth=2, color='blue')

plt.plot(1., 0., 'bo', markersize=4)
plt.plot(3., 0., 'bo', markersize=4)
plt.plot(5., 0., 'bo', markersize=4)

plt.axvline(2., linewidth=1, linestyle='dotted', color='black')
plt.axvline(4., linewidth=1, linestyle='dotted', color='black')
plt.axhline(0., linewidth=1, linestyle='dotted', color='black')
plt.axhline(0.5, linewidth=1, linestyle='dotted', color='black')
plt.axhline(-0.5, linewidth=1, linestyle='dotted', color='black')

plt.xlabel(r'$\xi/\pi$', fontsize=16)
plt.ylabel(r'$\hat{x}/\hat{W}$', fontsize=16)
plt.xlim(0., 6,)
plt.ylim(-1.5,1.5)

plt.savefig("Island.pdf")

#plt.show()
