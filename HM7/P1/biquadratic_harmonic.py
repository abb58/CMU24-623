#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt

#-------------------------------------------

plt.figure(1)
x=np.arange(-1.5, 1.5, 0.001)
U=x**4-(2*x**2)+1
x0=-1.0
U_A=(4*(x-x0)**2)
plt.plot(x, U, lw=2.0, label='$U(x)$')
plt.hold(True)
plt.plot(x, U_A, lw=2.0, label='$U(x)^{harm}$')
plt.ylim([0,2])
plt.xlabel('$x$',    fontsize=16)
plt.ylabel('$U(x)$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('P1-a.png',dpi=300)
plt.show()

#-------------------------------------------

plt.figure(2)
beta=np.arange(0.01, 100, 0.01)
k_htst=(1/(2*np.pi)) * np.sqrt(8) * np.exp(-beta)
plt.loglog(beta, k_htst, lw=2.0)
plt.hold(True)
beta_max=0.81
k_htst_max=(1/(2*np.pi)) * np.sqrt(8) * np.exp(-beta_max)
plt.loglog(beta_max, k_htst_max, 'ko', ms=8, label='limit')
plt.xlabel('$\\beta$',    fontsize=16)
plt.ylabel('$k^{h-TST}$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('P1-b.png',dpi=300)
plt.show()

#-------------------------------------------

