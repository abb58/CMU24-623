#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt


plt.figure(1)
beta=np.arange(0.01, 100, 0.01)
k_htst=(1/(2*np.pi)) * np.sqrt(8) * np.exp(-beta)
beta_mc=[0.01, 0.1, 1, 10, 100]

plt.loglog(beta, k_htst, lw=2.0)
plt.hold(True)



plt.xlabel('$\\beta$',    fontsize=16)
plt.ylabel('$k^{h-TST}$', fontsize=16)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('P1-b.png',dpi=300)
plt.show()



