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
t = data[:,1]
T = data[:,2]
P = data[:,3]
PE= data[:,4]
KE= data[:,5]
TE= data[:,6]
px= data[:,7]
py= data[:,8]
pz= data[:,9]

# Energy
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t, PE, '-b', lw=2.0, label='PE')
plt.hold(True)
plt.plot(t, KE, '-r', lw=2.0, label='KE')
plt.hold(True)
plt.plot(t, TE, '-g', lw=2.0, label='TE')
plt.hold(True)
plt.title('LJ Argon Liquid Simulation (Energy)', fontsize=18)
plt.xlabel('Time, t [units]', fontsize=16)
plt.ylabel('Energy', fontsize=16)
# Energy Drift
plt.subplot(2,1,2)
TE_cavg = moving_avg(data[10:,6])
drift = (data[10:,6] - TE_cavg)/256.0
plt.plot(data[10:,1], drift, '-k', lw=2.0, label='Energy Drift')
plt.hold(True)
drift_avg = moving_avg(drift)
print(drift_avg)
plt.plot(data[10:,1], drift_avg, '--r', lw=2.0, label='Avg_Energy Drift')
plt.xlabel('Time, t [units]', fontsize=16)
plt.ylabel('Energy Drift per atom', fontsize=16)
plt.legend(loc='best', fontsize=16)
plt.savefig('LJ-md-Ener.png',dpi=300)
plt.show()

# Temperature
fig, ax=plt.subplots()
t_sample=data[10:,1]
T_sample=data[10:,2]
T_cavg=moving_avg(T_sample)
ax.plot(t, T, '-b', lw=2.0, label='inst.T')
plt.hold(True)
ax.plot(t_sample, T_cavg, '-r', lw=2.0, label='Avg.T')
ax.ticklabel_format(useOffset=False)
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (Temperature)', fontsize=18)
plt.xlabel('Time, t [units]',fontsize=16)
plt.ylabel('Temperature [K]',fontsize=16)
plt.savefig('LJ-md-Temp.png',dpi=300)
plt.show()

# Pressure
plt.figure(3)
t_sample=data[10:,1]
P_sample=data[10:,3]
P_cavg=moving_avg(P_sample)
print(P_cavg)
plt.plot(t, P, '-b', lw=2.0, label='inst.P')
plt.hold(True)
plt.plot(t_sample, P_cavg, '-r', lw=2.0, label='Avg.P')
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (Pressure)', fontsize=18)
plt.xlabel('Time, t [units]',fontsize=16)
plt.ylabel('Pressure [MPa]',fontsize=16)
plt.savefig('LJ-md-Pressure.png',dpi=300)
plt.show()

# Momentum
plt.figure(4)
plt.plot(t, px, '-b', lw=2.0, label='Px')
plt.hold(True)
plt.plot(t, py, 'or', lw=2.0, label='Py')
plt.hold(True)
plt.plot(t, pz, '-g', lw=2.0, label='Pz')
plt.hold(True)
plt.legend(loc='best', fontsize=16)
plt.title('LJ Argon Liquid Simulation (Momentum)', fontsize=18)
plt.xlabel('Time, t [units]', fontsize=16)
plt.ylabel('Momentum',fontsize=16)
plt.savefig('LJ-md-Mom.png',dpi=300)
plt.show()
