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

sol = root( lambda T: CuAg_liq.G(T, (1,0) ) - CuAg_fcc.G(T, (1,0) ), 1000.)
print(sol.message)
T_Ag = sol.x[0]
print("melting point", T_Ag-T0)
plt.plot([1], [T_Ag], "o")

sol = root( lambda T: CuAg_liq.G(T, (0,1) ) - CuAg_fcc.G(T, (0,1) ), 1000.)
print(sol.message)
T_Cu = sol.x[0]
print("melting point", T_Cu-T0)
plt.plot([0], [T_Cu], "o")



from autograd import deriv, grad, hessian

def _grad(*a, **kw):
    f = grad(*a, **kw)
    return lambda T, y: np.array(f(T, np.array(y)))

def _hessian(*a, **kw):
    f = hessian(*a, **kw)
    return lambda T, y: f(T, np.array(y))

CuAg_fcc.dGdT = deriv(CuAg_fcc.G, argnum=0)
CuAg_liq.dGdT = deriv(CuAg_liq.G, argnum=0)

CuAg_fcc.dGdy = _grad(CuAg_fcc.G, argnum=1)
CuAg_liq.dGdy = _grad(CuAg_liq.G, argnum=1)

CuAg_fcc.d2GdTy = _grad(CuAg_fcc.dGdT, argnum=1)
CuAg_liq.d2GdTy = _grad(CuAg_liq.dGdT, argnum=1)

CuAg_fcc.d2Gdyy = _hessian(CuAg_fcc.G, argnum=1)
CuAg_liq.d2Gdyy = _hessian(CuAg_liq.G, argnum=1)


def MU(p, T, y):
    return p.G(T, (y, 1-y)) + p.dGdy(T, (y,1-y)) - np.array([y,1-y]) @ p.dGdy(T, (y,1-y))

T, yl, y1, y2 = 1000, 0.5, 0.05, 0.9


def func(p):
    print(".", end="", flush=True)
    (T, yl, y1, y2), mu = p[:4], p[4:]
    return np.concatenate([
        MU(CuAg_liq, T, yl) - mu,
        MU(CuAg_fcc, T, y1) - mu,
        MU(CuAg_fcc, T, y2) - mu,
    ])

sol = root(func, [T, yl, y1, y2, -55000, -45000], options=dict(factor=0.001))
print(sol.message)
T, yl, y1, y2 = sol.x[:4]
print("eutectic point", T-T0, yl, y1, y2)

plt.figure("diagram")
plt.plot([y1, yl, y2], [T, T, T], "o-")

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


def solve_line(M, v):
    n = len(v)
    A = np.ones((n+1,n+1))
    A[-1,-1] = 0
    A[:-1,:-1] = M

    b = np.zeros(n+1)
    b[:-1] = v
    
    return np.linalg.solve(A, b)[:-1]

def func(T, p, phA, phB):
    print(".", end="", flush=True)
    y1A, y2A, y1B, y2B = p 
    dmu = np.linalg.solve( #fill )

    dyA = solve_line( #fill)
    dyB = solve_line( #fill)

    return np.concatenate([ dyA, dyB ])

from scipy.integrate import solve_ivp

sol = solve_ivp(func, (T, 500), [#fill], args=(#fill) , max_step=5)
print()
print(sol.message)
plt.figure("diagram")
plt.plot(sol.y[0], sol.t, "-")
plt.plot(sol.y[2], sol.t, "-")

sol = solve_ivp(func, (T, T_Cu), [#fill], args=(#fill) , max_step=5)
print()
print(sol.message)
plt.figure("diagram")
plt.plot(sol.y[0], sol.t, "-")
plt.plot(sol.y[2], sol.t, "-")

sol = solve_ivp(func, (T, T_Ag), [#fill], args=(#fill) , max_step=5)
print()
print(sol.message)
plt.figure("diagram")
plt.plot(sol.y[0], sol.t, "-")
plt.plot(sol.y[2], sol.t, "-")


plt.show()
