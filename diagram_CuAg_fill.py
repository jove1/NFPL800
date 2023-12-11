#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

y = np.linspace(0, 1, 250)

from data_CuAg import CuAg_liq, CuAg_fcc
CuAg_liq.y = y, 1-y
CuAg_fcc.y = y, 1-y

phases = [CuAg_liq, CuAg_fcc]

plt.figure("diagram")

N = sum( p.X(p.y)[0].size for p in phases )
points = np.empty( (N, 2) )
phase_id = np.empty(N, dtype=int)

for T in range(500, 1500, 15):
    print(".", end="", flush=True)

    o = 0
    for i, p in enumerate(phases):
        X = p.X(p.y)[0]
        points[o:o+X.size, 0] = X
        points[o:o+X.size, 1] = p.G(T, p.y)
        phase_id[o:o+X.size] = i
        o += X.size

    hull = ConvexHull(points)
    mask = hull.equations[:,1] < 0 
    for a, b in hull.simplices[mask]:
        if phase_id[a] != phase_id[b] or abs(points[a,0] - points[b,0]) > 0.05:
            plt.plot([ points[a, 0], points[b, 0] ], [ T, T ], "k-")
print()
plt.xlim(0,1)


from scipy.optimize import root

T0 = 273.15

sol = root( lambda T: # fill in here )
print(sol.message)
T = sol.x[0]
print("Ag melting point", T-T0)
plt.plot([1], [T], "o")

sol = # fill this as well for melting point of Cu
print(sol.message)
T = sol.x[0]
print("Cu melting point", T-T0)
plt.plot([0], [T], "o")



from data_CuAg import GHSERAG, GHSERCU, R
liq_dGdy = lambda T, y: (
    # calculate derivative of CuAg_liq.G (data_CuAg.py) with respect to y  by hand
)


from autograd import grad
CuAg_fcc.dGdy = grad(CuAg_fcc.G, argnum=1)
CuAg_liq.dGdy = grad(CuAg_liq.G, argnum=1)

T = 1000

plt.figure("test derivative")
y = np.linspace(0, 1, 250)[1:-1]
plt.xlim(0,1)
plt.plot(y, [CuAg_liq.dGdy(T, (_y, 1-_y) )[0] for _y in y] )
plt.plot(y, [CuAg_liq.dGdy(T, (_y, 1-_y) )[1] for _y in y] )

plt.plot(y, liq_dGdy(T, (y, 1-y) )[0], ".")
plt.plot(y, liq_dGdy(T, (y, 1-y) )[1], ".")


def MU(p, T, y):
    # how to calculate chemical potential of phase p at T and y?
    return # fill this ...

T, yl, y1, y2 = 1000, 0.5, 0.05, 0.9

# 
#def func(p):
#    print(".", end="", flush=True)
#    (T, yl, y1, y2), mu = p[:4], p[4:]
#    return np.concatenate([
#    # then figure what to put here to get 3-phase equilibrium ... 
#   ])
#
#sol = root(func, [T, yl, y1, y2, -55000, -45000], options=dict(factor=0.001))
#print(sol.message)
#T, yl, y1, y2 = sol.x[:4]
#print("eutectic point", T-T0, yl, y1, y2)

#plt.figure("diagram")
#plt.plot([y1, yl, y2], [T, T, T], "o-")
# 

plt.figure("Gibbs")
plt.xlim(0,1)
y = np.linspace(0, 1, 250)
plt.plot(y, CuAg_liq.G(T, (y, 1-y)) )
plt.plot(y, CuAg_fcc.G(T, (y, 1-y)) )

plt.plot([1,0], MU(CuAg_liq, T, yl), "o--")
plt.plot([1,0], MU(CuAg_fcc, T, y1), "o--")
plt.plot([1,0], MU(CuAg_fcc, T, y2), "o--")

plt.plot([yl, y1, y2], [
    CuAg_liq.G(T, (yl, 1-yl)),
    CuAg_fcc.G(T, (y1, 1-y1)),
    CuAg_fcc.G(T, (y2, 1-y2)),
], "o")

plt.show()
