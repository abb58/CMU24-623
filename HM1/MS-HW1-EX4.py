import matplotlib
import numpy as np
import matplotlib.pylab as plt

#-----------------------------------------------------
def plot_lj1(Rnn):
    return( 4*((1/Rnn)**12 - (1/Rnn)**6) )

def plot_lj2(Rnn,A12,A6):
    return( 2*(A12*(1/Rnn)**12 - A6*(1/Rnn)**6) )
#-----------------------------------------------------

grid = 100
r    = np.linspace(0.90, 1.75, grid)
A12  = 12.13
A6   = 14.45

U_r   = plot_lj1(r)
U_rnn = plot_lj2(r,A12,A6)

plt.plot(r, np.zeros(len(r)), '-k')
plt.hold(True)
plt.plot(r, U_r, label='Eq:1')
plt.hold(True)
plt.plot(r, U_rnn, label='Eq:3')
plt.hold(True)
plt.plot( 1, 0,'-bo', label='Eq:1=0')
plt.hold(True)
plt.plot( 0.97125, 0,'-go', label='Eq:3=0')
plt.hold(True)
plt.plot( 1.12246, min(U_r),'-bd', label='min(Eq:1)')
plt.hold(True)
plt.plot( 1.09019, min(U_rnn),'-gd', label='min(Eq:3)')

plt.xlabel("$R_{nn}$")
plt.ylabel("$U(r_{nn}$)")
plt.legend(loc='best')
plt.savefig('LJ-Ex4.png',dpi=300)
plt.show()
