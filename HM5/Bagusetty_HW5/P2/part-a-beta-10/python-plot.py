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

start = 10

data1=np.loadtxt('LDmj_sim_dx-0.05.txt')
t1 = data1[start:,0]
U1 = data1[start:,1]
x1 = data1[start:,2]
x12= data1[start:,3]
data2=np.loadtxt('LDmj_sim_dx-0.5.txt')
t2 = data2[start:,0]
U2 = data2[start:,1]
x2 = data2[start:,2]
x22= data2[start:,3]
data3=np.loadtxt('LDmj_sim_dx-5.txt')
t3 = data3[start:,0]
U3 = data3[start:,1]
x3 = data3[start:,2]
x32= data3[start:,3]

# Potential
U1_avg = moving_avg(data1[start:,1])
U2_avg = moving_avg(data2[start:,1])
U3_avg = moving_avg(data3[start:,1])

plt.figure(1)
plt.plot(t1, U1_avg, '-r', lw=2.0, label='$dx$=0.05')
plt.hold(True)
plt.plot(t2, U2_avg, '-g', lw=2.0, label='$dx$=0.5')
plt.hold(True)
plt.plot(t3, U3_avg, '-b', lw=2.0, label='$dx$=5')
plt.hold(True)
plt.title('Metropolis MC Simulation (U), $ \\beta=10$', fontsize=18)
plt.xlabel('Moves',    fontsize=16)
plt.ylabel('$\\langle U \\rangle$',        fontsize=16)
plt.legend(loc='best', fontsize=16)
plt.savefig('Oscillator2a_U.png',dpi=300)
plt.show()

# Position
x1_avg = moving_avg(data1[start:,2])
x2_avg = moving_avg(data2[start:,2])
x3_avg = moving_avg(data3[start:,2])

x21_avg = moving_avg(data1[start:,3])
x22_avg = moving_avg(data2[start:,3])
x23_avg = moving_avg(data3[start:,3])

fig, ax=plt.subplots(2, sharex=True)
ax[0].set_title('Metropolis MC Simulation ($x, x^2$), $\\beta=10$', fontsize=18)

ax[0].plot(t1, x1_avg, '-r', lw=2.0, label='$dx$=0.05')
plt.hold(True)
ax[0].plot(t2, x2_avg, '-g', lw=2.0, label='$dx$=0.5')
plt.hold(True)
ax[0].plot(t3, x3_avg, '-b', lw=2.0, label='$dx$=5')
plt.hold(True)
ax[0].set_ylabel('$\\langle x \\rangle$',     fontsize=16)
plt.legend(loc='best', fontsize=16)


ax[1].plot(t1, x21_avg, '-r', lw=2.0, label='$dx$=0.05')
plt.hold(True)
ax[1].plot(t2, x22_avg, '-g', lw=2.0, label='$dx$=0.5')
plt.hold(True)
ax[1].plot(t3, x23_avg, '-b', lw=2.0, label='$dx$=5')
plt.hold(True)
plt.legend(loc='best', fontsize=16)
plt.hold(True)
plt.xlabel('Moves', fontsize=16)
plt.ylabel('$\\langle x^2 \\rangle$',     fontsize=16)
plt.legend(loc='best', fontsize=16)

plt.savefig('Oscillator2a-x.png',dpi=300)
plt.show()
