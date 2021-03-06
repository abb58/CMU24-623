#/usr/bin/python


from scipy.integrate import quad
from pylab import *
beta=[0.1, 1.0, 5.0, 10.0]


for i in beta:
    def potential1(x):
        return (0.5*x**2)

    def potential2(x):
        return (x**4 - (2*x**2) + 1)

    def integrand1(x):
        return (x * exp(-i*potential1(x)))

    def integrand2(x):
        return exp(-i*potential1(x))
    
    ans1, err1 = quad(integrand1, -inf, inf)
    ans2, err2 = quad(integrand2, -inf, inf)
    print ans1/ans2
