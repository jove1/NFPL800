#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigureSlope

np.seterr(all='print')

R = 8.314

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

class Phase:
    def __init__(self, Ga, Gb, L0, L1):
        self.Ga = Ga
        self.Gb = Gb
        self.L0 = L0
        self.L1 = L1

    def G(self, x):
        xa, xb = x
        return self.Ga*xa + self.Gb*xb + ( self.L0 + self.L1*(xa-xb) )*xa*xb + RT*( xlnx(xa) + xlnx(xb) )

    def dGdx(self, x):
        xa, xb = x
#PROBLEM
        # Write gradient of self.G as np.array([ dG/dx[0], dG/dx[1] ])
        return np.array([0, 0]) 
#END

I = np.eye(2)
RT = R * 300
Ga = 500
Gb = 2000
L0 = -4000
#L0 = 6000
L1 = 2000

phase = Phase(Ga, Gb, L0, L1)

def calc_root(slope):
    from scipy import optimize
#PROBLEM
    # Write function that has root at a given slope of phase.G,
    # don't forget to include constraint xa + xb = 1
    def func(x):
       return x - 0.5
#END

    sol = optimize.root(func, 
        x0=[0.5, 0.5],
        options={'factor': 1e-3}
    )
    print(sol.x, sol.message)
    x = sol.x

    g = phase.G(x)
#PROBLEM
    # Calculate chemical potential at correct offset for phase.G
    mu = np.array([0, slope])
#END
    return [0, x[1], 1], [mu[0], g, mu[1]]

with InteractiveFigureSlope("phase.G([xa,xb]), root", calc_root) as ax:
    x = np.linspace(0, 1, 150)
    ax.plot(x, phase.G([1-x, x]))

def calc_minimize(slope):
    from scipy import optimize
#PROBLEM
    # Find given slope of phase.G as a minimum of a function,
    # don't forget to include correct constraint xa + xb = 1
    sol = optimize.minimize(lambda x: (x[0]-0.5)**2,
        x0=[0.5, 0.5],
        bounds=[(0, 1)]*2,
        constraints = [ dict(type="eq", fun=lambda x: x[1] - .5) ]
    )
#END
    print(sol.x, sol.message)
    x = sol.x

    g = phase.G(x)
#PROBLEM
    # Calculate chemical potential at correct offset for phase.G
    mu = np.array([0, slope])
#END
    return [0, x[1], 1], [mu[0], g, mu[1]]

with InteractiveFigureSlope("phase.G([xa,xb]), minimize", calc_minimize) as ax:
    x = np.linspace(0, 1, 150)
    ax.plot(x, phase.G([1-x, x]))

plt.show()
