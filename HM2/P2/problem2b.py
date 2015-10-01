#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt

Vo = 1.414;
t = np.linspace(0,20,160)
x = Vo*np.sin(t)
v = Vo*np.cos(t)

# Problem 1
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(t,x,'-k',linewidth=2.0,label='position')
plt.hold(True)
plt.plot(t,v,'--k',linewidth=2.0, label='velocity')
plt.xlabel('t')
plt.ylabel('x(t) / v(t)')
plt.title('Case (i) Linear Spring')
plt.legend(loc='best')
plt.show()

#----------- Problem 2 ------------
x1c,v1c=[],[]
dt = 0.125  # time step
N = 20/dt # No of steps
t = np.linspace(0,20,N) # time  array

# Velocity Verlet Scheme
for i in range(len(t)):
    tempx = Vo*( np.sin(t[i]) + np.cos(t[i])*dt - 0.5*dt**2*np.sin(t[i]) )
    x1c.append(tempx)

    tempv = Vo*( np.cos(t[i]) - 0.5*dt*np.sin(t[i]) - 0.5*dt*np.sin(t[i])*np.cos(dt) - 0.5*dt*np.cos(t[i])*np.sin(dt) )
    v1c.append(tempv)

# Conservation of Energy
Ef = [0.5*i*i for i in x1c]
kf = [0.5*j*j for j in v1c]
Tf = [a + b for a, b in zip(Ef, kf)]
print(len(Tf))
print(len(t))
Ti = np.ones(len(Tf))
plt.plot(t, Tf, '-k', label='Energy-final')
plt.hold(True)
plt.plot(t, Ti, 'ok', label='Energy-initial')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Case (i) Linear Spring (VV scheme vs Analytical)')
plt.legend(loc='best')
plt.show()

    
fig1c_pos = plt.figure()
ax = fig1c_pos.add_subplot(111)
plt.plot(t,x1c,'-k',linewidth=2.0,label='position-VV')
plt.hold(True)
plt.plot(t,x,'ok',linewidth=1.0, label='position-(1b)')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Case (i) Linear Spring (VV scheme vs Analytical)')
plt.legend(loc='best')
plt.show()

fig1c_vel = plt.figure()
ax = fig1c_vel.add_subplot(111)
plt.plot(t,v1c,'-k',linewidth=2.0,label='Velocity-VV')
plt.hold(True)
plt.plot(t,v,'dk',linewidth=1.0, label='Velocity-(1b)')
plt.xlabel('t')
plt.ylabel('v(t)')
plt.title('Case (i) Linear Spring (VV scheme vs Analytical)')
plt.legend(loc='best')
plt.show()

fig1c_posvel = plt.figure()
ax = fig1c_posvel.add_subplot(111)
plt.plot(x1c,v1c,'-k',linewidth=2.0)
plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Case (i) Linear Spring')
plt.legend(loc='best')
plt.show()


x = np.linspace(-4,4,160)
Usa = 0.5*x*x
Usb = [0.5*i*i for i in x1c]
fig1c_pot = plt.figure()
ax = fig1c_pot.add_subplot(111)
plt.plot(x,Usa,'-k',linewidth=1.0)
plt.hold(True)
plt.plot(x1c,Usb,'xk')
plt.xlabel('x(t)')
plt.ylabel('U(t)')
plt.title('Case (i) Linear Spring')
plt.legend(loc='best')
plt.show()

#----------- Problem 1d ------------
x1c,v1c=[],[]
dt = 0.2  # time step
N = 20/dt # No of steps
t = np.linspace(0,20,N) # time  array

# Velocity Verlet Scheme
for i in range(len(t)):
    tempx = Vo*( np.sin(t[i]) + np.cos(t[i])*dt - 0.5*dt**2*np.sin(t[i]) )
    x1c.append(tempx)

    tempv = Vo*( np.cos(t[i]) - 0.5*dt*np.sin(t[i]) - 0.5*dt*np.sin(t[i])*np.cos(dt) - 0.5*dt*np.cos(t[i])*np.sin(dt) )
    v1c.append(tempv)

    
fig1c_pos = plt.figure()
ax = fig1c_pos.add_subplot(111)
plt.plot(t,x1c,'-k',linewidth=2.0,label='position-VV')
plt.hold(True)
plt.plot(t,x,'ok',linewidth=1.0, label='position-(1b)')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Case (i) Linear Spring (VV scheme vs Analytical)')
plt.legend(loc='best')
plt.show()

fig1c_vel = plt.figure()
ax = fig1c_vel.add_subplot(111)
plt.plot(t,v1c,'-k',linewidth=2.0,label='Velocity-VV')
plt.hold(True)
plt.plot(t,v,'dk',linewidth=1.0, label='Velocity-(1b)')
plt.xlabel('t')
plt.ylabel('v(t)')
plt.title('Case (i) Linear Spring (VV scheme vs Analytical)')
plt.legend(loc='best')
plt.show()

fig1c_posvel = plt.figure()
ax = fig1c_posvel.add_subplot(111)
plt.plot(x1c,v1c,'-k',linewidth=2.0)
plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Case (i) Linear Spring')
plt.legend(loc='best')
plt.show()


x = np.linspace(-4,4,100)
Usa = 0.5*x*x
Usb = [0.5*i*i for i in x1c]
fig1c_pot = plt.figure()
ax = fig1c_pot.add_subplot(111)
plt.plot(x,Usa,'-k',linewidth=1.0)
plt.hold(True)
plt.plot(x1c,Usb,'xk')
plt.xlabel('x(t)')
plt.ylabel('U(t)')
plt.title('Case (i) Linear Spring')
plt.legend(loc='best')
plt.show()

