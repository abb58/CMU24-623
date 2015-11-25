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
PE= data[:,2]

# Energy
plt.figure(1)
t_sample=data[10:,0]
U_sample=data[10:,2]
U_cavg=moving_avg(P_sample)
plt.plot(t, PE, '-b', lw=2.0, label='PE')
plt.hold(True)
plt.plot(t_sample, U_cavg, '-r', lw=2.0, label='Avg.P')
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (Energy)', fontsize=18)
plt.xlabel('Time, t [units]', fontsize=16)
plt.ylabel('U', fontsize=16)
plt.savefig('LJ-md-Ener.png',dpi=300)
plt.show()

# Pressure
plt.figure(2)
t_sample=data[10:,0]
P_sample=data[10:,1]
P_cavg=moving_avg(P_sample)
print(P_cavg)
plt.plot(t, P, '-b', lw=2.0, label='inst.P')
plt.hold(True)
plt.plot(t_sample, P_cavg, '-r', lw=2.0, label='Avg.P')
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (Pressure)', fontsize=18)
plt.xlabel('Time, t [units]',fontsize=16)
plt.ylabel('dim.less Pressure',fontsize=16)
plt.savefig('LJ-md-Pressure.png',dpi=300)
plt.show()
