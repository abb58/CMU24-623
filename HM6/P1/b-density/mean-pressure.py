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

data =np.loadtxt('LDmj_sim_7.712.ener')

P_sample=data[10:,1]
t_sample=data[10:,0]
P_cavg=moving_avg(P_sample)
print(P_cavg)
