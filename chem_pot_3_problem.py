#!/usr/bin/env python3
# Gibbs-energy models and chemical potentials for phase-equilibrium exercises.

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigureSlope

np.seterr(all='print')

R = 8.314  # Gas constant, in J/(mol K).

# Ideal-mixing contribution, with the limiting value x*ln(x) = 0 at x = 0.
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

# L0 and L1 set the symmetric and composition-asymmetric interaction contributions.
L0 = -4000
#L0 = 6000
L1 = 2000
RT = R * 300

# Binary molar Gibbs energy, including a higher-order interaction term (L1).
def G(x):
    return L0*x*(1-x) + L1*x*(1-x)*(1-2*x) + RT*( xlnx(x) + xlnx(1-x) )

# The derivative gives the local slope of G(x).
def dGdx(x):
#PROBLEM
    return 0 # Write derivative of G
#END

x = np.linspace(0, 1, 150)

plt.figure()
plt.title("G(x), root")
plt.xlim(0, 1)
plt.plot(x, G(x))

# The callback receives a target tangent slope and returns coordinates for plotting.
@InteractiveFigureSlope
def calc_root(slope):
    from scipy import optimize
#PROBLEM
    # Find x for a given slope of G as a root of an equation
    sol = optimize.root(lambda x: x-0.5, 
        x0=[0.5],
        #options={'factor': 1e-3}
    )
#END
    print(sol.x, sol.message)

    x = sol.x[0]
    g = G(x)
    return [0, x, 1], [g - x*slope, g, g + (1-x)*slope]

plt.figure()
plt.title("G(x), minimize")
plt.xlim(0, 1)
plt.plot(x, G(x))

@InteractiveFigureSlope
def calc_minimize(slope):
    from scipy import optimize
#PROBLEM
    # Find x for a given slope of G as a minimum of a function
    sol = optimize.minimize(lambda x: (x-0.5)**2, 
        x0=[0.5],
        method='Nelder-Mead',
        #bounds=[(0, 1)]
    )
#END
    print(sol.x, sol.message)

    x = sol.x[0]
    g = G(x)
    return [0, x, 1], [g - x*slope, g, g + (1-x)*slope]

plt.show()
