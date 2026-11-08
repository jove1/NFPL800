#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigurePos

np.seterr(all='print')

R = 8.314

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

L = -4000
#L = 6000
RT = R * 300

def G(x):
    return L*x*(1-x) + RT*( xlnx(x) + xlnx(1-x) )

def dGdx(x):
#PROBLEM
    return 0 # Write derivative of G
#END

x = np.linspace(0, 1, 150)

plt.figure()
plt.title("G(x)")
plt.xlim(0, 1)
plt.plot(x, G(x))

@InteractiveFigurePos
def calc(x):
    g = G(x)
    dg = dGdx(x)
#PROBLEM 
    return [0, x, 1], [
        g, # Write value of tanget at 0 (chemical potential)
        g, 
        g  # Write value of tanget at 1 (chemical potential)
    ]
#END

plt.show()
