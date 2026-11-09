#!/usr/bin/env python3
# Gibbs-energy models and chemical potentials for phase-equilibrium exercises.

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigureSlope

np.seterr(all='print')

R = 8.314  # Gas constant, in J/(mol K).

# Ideal-mixing contribution, with x*ln(x) defined as 0 in the limit x -> 0.
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

# Binary phase with end-member energies and composition-dependent interaction terms.
class Phase:
    def __init__(self, Ga, Gb, L0, L1):
        self.Ga = Ga
        self.Gb = Gb
        self.L0 = L0
        self.L1 = L1

    # Molar Gibbs energy for the two component fractions [xa, xb].
    def G(self, x):
        xa, xb = x
        return self.Ga*xa + self.Gb*xb + ( self.L0 + self.L1*(xa-xb) )*xa*xb + RT*( xlnx(xa) + xlnx(xb) )

    # Gradient of G with respect to both independent composition variables.
    def dGdx(self, x):
        xa, xb = x
#PROBLEM
        # Write gradient of self.G as np.array([ dG/dx[0], dG/dx[1] ])
        return np.array([0, 0]) 
#END

I = np.eye(2) # Identity matrix
RT = R * 300
Ga = 500
Gb = 2000
L0 = -4000
#L0 = 6000
L1 = 2000

phase = Phase(Ga, Gb, L0, L1)

# Evaluate the energy along xa + xb = 1 for plotting.
x = np.linspace(0, 1, 150)

plt.figure()
plt.title("phase.G([xa,xb]), root")
plt.xlim(0, 1)
plt.plot(x, phase.G([1-x, x]))

# The callback receives a requested slope and returns three points on the tangent.
# The interactive module only consumes the returned x- and y-coordinate sequences.
@InteractiveFigureSlope
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


x = np.linspace(0, 1, 150)

# The minimization approach
plt.figure()
plt.title("phase.G([xa,xb]), minimize")
plt.xlim(0, 1)
plt.plot(x, phase.G([1-x, x]))

@InteractiveFigureSlope
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

plt.show()
