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

data=np.loadtxt('LDmj_sim.ener')
t = data[:,0]
P = data[:,1]
U = data[:,2]

# Energy
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t, U, '-b', lw=2.0, label='PE')
plt.title('LJ Argon Liquid Simulation (Energy)', fontsize=18)
plt.xlabel('Moves', fontsize=16)
plt.ylabel('U', fontsize=16)

# Energy Drift
plt.subplot(2,1,2)
t_drift = data[1000:,0]
U_sample= data[1000:,2]
U_avg   = np.mean(U_sample)
print('Avg of U : {0}'.format(U_avg))
U_drift = (U_sample - U_avg)**2
print('Value of drift : {0}'.format(np.mean(U_drift)))

plt.plot(t_drift, U_drift, '-r', lw=2.0, label='Energy Drift')
plt.xlabel('Moves', fontsize=16)
plt.ylabel('(E-$\\langle E \\rangle)^2$', fontsize=16)
plt.legend(loc='best', fontsize=16)
plt.savefig('LJ-md-Ener.png',dpi=300)
plt.show()

# # Pressure
# plt.figure(2)
# t_sample=data[10:,0]
# P_sample=data[10:,1]
# P_cavg=moving_avg(P_sample)
# print(P_cavg)
# plt.plot(t, P, '-b', lw=2.0, label='inst.P')
# plt.hold(True)
# plt.plot(t_sample, P_cavg, '-r', lw=2.0, label='Avg.P')
# plt.legend(loc='best', fontsize=16)
# plt.title('LJ Argon Liquid Simulation (Pressure)', fontsize=18)
# plt.xlabel('Time, t [units]',fontsize=16)
# plt.ylabel('Pressure [MPa]',fontsize=16)
# plt.savefig('LJ-md-Pressure.png',dpi=300)
# plt.show()
