#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigurePos

np.seterr(all='print')

R = 8.314

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

class Phase:
    def __init__(self, Ga, Gb, L):
        self.Ga = Ga
        self.Gb = Gb
        self.L = L

    def G(self, x):
        xa, xb = x
        return self.Ga*xa + self.Gb*xb + self.L*xa*xb + RT*( xlnx(xa) + xlnx(xb) )

    def dGdx(self, x):
        xa, xb = x
#PROBLEM
        # Write gradient of self.G as np.array([ dG/dxa, dG/dxb ])
        return np.array([0, 0]) 
#END

I = np.eye(2)
RT = R * 300
Ga = 500
Gb = 2000
L = -4000
#L = 6000

phase = Phase(Ga, Gb, L)

def calc(x1):
    x = np.array([1-x1, x1])
    g = phase.G(x)
    dg = phase.dGdx(x)
#PROBLEM
    return [0, x1, 1], [
        g, # Write value of chemical potential
        g, 
        g  # Write value of chemical potential
    ]
#END

with InteractiveFigurePos("phase.G([xa, xb])", calc) as ax:
    x = np.linspace(0, 1, 150)
    ax.plot(x, phase.G([1-x, x]))

plt.show()
