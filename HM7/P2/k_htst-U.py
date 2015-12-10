#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt


plt.figure(1)
beta=np.arange(0.01, 100, 0.01)
k_htst=(1/(2*np.pi)) * np.sqrt(8) * np.exp(-beta)


plt.loglog(beta, k_htst, lw=2.0)
plt.hold(True)
plt.loglog(beta_mc, k_mc_htst, 'ko', ms=10)



plt.xlabel('$\\beta$',    fontsize=16)
plt.ylabel('$k^{h-TST}$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('P2-maxstepsize.png',dpi=300)
plt.show()



