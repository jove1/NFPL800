#!/usr/bin/env python3
# Gibbs-energy models and chemical potentials for phase-equilibrium exercises.

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigurePos

np.seterr(all='print')

R = 8.314  # Gas constant, in J/(mol K).

# The ideal-mixing contribution is defined by its limit, x*ln(x) -> 0 at x -> 0.
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

# Binary phase with end-member energies Ga and Gb and interaction parameter L.
class Phase:
    def __init__(self, Ga, Gb, L):
        self.Ga = Ga
        self.Gb = Gb
        self.L = L

    # Molar Gibbs energy as a function of the two site fractions [xa, xb].
    def G(self, x):
        xa, xb = x
        return self.Ga*xa + self.Gb*xb + self.L*xa*xb + RT*( xlnx(xa) + xlnx(xb) )

    # Gradient with respect to [xa, xb]; both components are needed for the tangent.
    def dGdx(self, x):
        xa, xb = x
#PROBLEM
        # Write gradient of self.G as np.array([ dG/dxa, dG/dxb ])
        return np.array([0, 0]) 
#END

I = np.eye(2) # Identity matrix
RT = R * 300  # Thermal energy at 300 K.
Ga = 500
Gb = 2000
L = -4000
#L = 6000

phase = Phase(Ga, Gb, L)

# Plot the phase energy along the binary-composition line xa + xb = 1.
x = np.linspace(0, 1, 150)

plt.figure()
plt.title("phase.G([xa, xb])")
plt.xlim(0, 1)
plt.plot(x, phase.G([1-x, x]))

# calc(x1) returns the x- and y-coordinates of the common-tangent construction.
@InteractiveFigurePos
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


plt.show()
