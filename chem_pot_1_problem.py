#!/usr/bin/env python3
# Gibbs-energy models and chemical potentials for phase-equilibrium exercises.

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigurePos

np.seterr(all='print')

R = 8.314  # Gas constant, in J/(mol K).

# The ideal-mixing term x*ln(x) has the limiting value 0 at x = 0.
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

# L controls the energetic preference for mixing; its sign and magnitude affect G(x).
L = -4000
#L = 6000
RT = R * 300

# Molar Gibbs energy of a symmetric binary regular solution.
def G(x):
    return L*x*(1-x) + RT*( xlnx(x) + xlnx(1-x) )

# Composition derivative: its value is the slope of the G(x) curve.
def dGdx(x):
#PROBLEM
    return 0 # Write derivative of G
#END

# Evaluate the free-energy curve over the full binary-composition range.
x = np.linspace(0, 1, 150)

plt.figure()
plt.title("G(x)")
plt.xlim(0, 1)
plt.plot(x, G(x))

# calc(x) returns plot coordinates and the corresponding tangent-line values as two arrays/lists.
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
