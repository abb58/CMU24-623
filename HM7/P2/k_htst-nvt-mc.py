#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt

plt.figure(1)
beta=np.arange(0.01, 100, 0.01)
k_htst=(1/(2*np.pi)) * np.sqrt(8) * np.exp(-beta)

beta_mc=[0.01, 0.1, 1, 10, 100]
#k_mc_htst=[2.43714, 1.51812, 0.633919, 0.200536, 0.0634151]
#k_mc_htst=[0.240044, 0.149243, 0.0610422, 0.0193045, 0.00610461]
k_mc_htst=[0.392001, 0.376149, 0.07155430, 1.26157E-05, 0]

plt.loglog(beta, k_htst, lw=2.0)
plt.hold(True)
plt.loglog(beta_mc, k_mc_htst, 'ro', ms=10)


plt.xlabel('$\\beta$',    fontsize=16)
plt.ylabel('$k^{h-TST}$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('P1-b.png',dpi=300)
plt.show()



