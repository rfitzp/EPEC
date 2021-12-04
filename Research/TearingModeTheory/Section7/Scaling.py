import matplotlib.pyplot as plt
import numpy as np

def rho1(nu):
    return nu

r1 = np.arange(0., 2., 1.e-3)

plt.hlines(y=1., xmin=0, xmax=2., linewidth=2, color='blue')
plt.plot(r1, rho1(r1), linewidth=2, color='blue')

plt.text(1., 0.5, r'${\rm VR}-\varphi$', fontsize=20, color='black')
plt.text(1.65, 1.5, r'${\rm VR}-\theta$', fontsize=20, color='black')
plt.text(0.1, 0.5, r'${\rm HR}-\varphi$', fontsize=20, color='black')
plt.text(0.5, 1.5, r'${\rm HR}-\theta$', fontsize=20, color='black')

plt.text(0.35, 1.04, r'$\hat{\rho}_\ast=\beta_\ast$', fontsize=20, color='red')
plt.text(0.3, 0.2, r'$\hat{\rho}_\ast=\beta_\ast^{-1}\hat{\nu}_\ast$', fontsize=20, color='red')

plt.xlabel(r'$\hat{\nu}_\ast/\beta_\ast^{\,2}$', fontsize=16)
plt.ylabel(r'$\hat{\rho}_\ast/\beta_\ast$', fontsize=16)

plt.savefig("Scaling.pdf")

#plt.show()
