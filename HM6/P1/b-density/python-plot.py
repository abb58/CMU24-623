#!/usr/bin/python

import matplotlib.pylab as plt
import numpy as np

def moving_avg(x):
    y = np.zeros(len(x))
    current_sum=0.0
    for i in range(len(x)):
        current_sum+= x[i]
        y[i] = current_sum/(i+1.0)
    return y

data =np.loadtxt('LDmj_sim_7.5070.ener')

# Energy
plt.figure(1)
t_sample=data[100:,0]; U_sample=data[100:,2];  U_cavg=moving_avg(U_sample)

plt.plot(t_sample, U_cavg, '-k', lw=2.0, label='MD-NVT')
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (U)', fontsize=18)
plt.xlabel('Moves', fontsize=16)
plt.ylabel('$\\langle U \\rangle$', fontsize=16)
plt.savefig('LJ-md-Ener.png',dpi=300)
plt.show()

# Pressure
plt.figure(2)
P_sample=data[100:,1]
P_cavg=moving_avg(P_sample)
print(P_cavg)
plt.plot(t_sample, P_sample, '-r', lw=2.0, label='inst. P')
plt.hold(True)
plt.plot(t_sample, P_cavg, '-b', lw=2.0, label='avg. P')
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (Pressure)', fontsize=18)
plt.xlabel('Moves',fontsize=16)
plt.ylabel('$\\langle P \\rangle$', fontsize=16)
plt.savefig('LJ-md-Pressure.png',dpi=300)
plt.show()
