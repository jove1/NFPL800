#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

R = 8.314
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

class Phase:
    pass

y = np.linspace(0, 1, 250)

a = Phase()
a.G = lambda T, y: y[0]*(-7000+120*T) + y[1]*(-6000+130*T) + R*T*( xlnx(y[0]) + xlnx(y[1]) ) + y[0]*y[1]*(35000-10*T -4000*(y[0]-y[1]))
a.X = lambda y: y
a.y = y, 1-y

b = Phase()
b.G = lambda T, y: y[0]*(-5000+120*T) + y[1]*(-8000+130*T) + R*T*( xlnx(y[0]) + xlnx(y[1]) ) + y[0]*y[1]*(30000-10*T)
b.X = lambda y: y
b.y = y, 1-y

c = Phase()
c.G = lambda T, y: y[0]*(+4000+110*T) + y[1]*(+5000+120*T) + R*T*( xlnx(y[0]) + xlnx(y[1]) ) + y[0]*y[1]*(20000-5*T)
c.X = lambda y: y
c.y = y, 1-y

d = Phase()
d.G = lambda T, y: 54000+ 119*(T-500)
d.X = lambda y: np.array([.5, .5])
d.y = None


from data_CuAg import CuAg_liq, CuAg_fcc
#CuAg_liq.y = #????
#CuAg_fcc.y = #????

#from data_MgNi import MgNi_liq, MgNi_fcc, MgNi_hcp, MgNi2, Mg2Ni
# ???

phases = [a, b, c, d]
#phases = [CuAg_liq, CuAg_fcc]
#phases = [MgNi_liq, MgNi_fcc, MgNi_hcp, MgNi2, Mg2Ni]


N = sum( p.X(p.y)[0].size for p in phases )
points = np.empty( (N, 2) )
phase_id = np.empty(N, dtype=int)

fig, ax = plt.subplots(1, 2, figsize=(14,6))
ax[1].grid()
ax[1].set_ylim(500, 1500)
ax[1].set_xlim(0, 1)

for T in range(500, 1500, 10):
    print(".", end="", flush=True)

    ax[0].cla()
    ax[0].grid()
    ax[0].set_xlim(0, 1)
    for p in phases:
        #ax[0].plot(????, ?????, ".")

    o = 0
    for i, p in enumerate(phases):
        X = #????
        points[o:o+X.size, 0] = X
        points[o:o+X.size, 1] = # ????
        phase_id[o:o+X.size] = i
        o += X.size

    hull = ConvexHull(points)
    mask = hull.equations[:,1] < 0 
    for a, b in hull.simplices[mask]:
        if phase_id[a] != phase_id[b]: #?????
            ax[1].plot([ points[a, 0], points[b, 0] ], [ T, T ], "k-")
            ax[1].plot([ points[a, 0] ], [ T ], f"C{phase_id[a]}.")
            ax[1].plot([ points[b, 0] ], [ T ], f"C{phase_id[b]}.")
            ax[0].plot([ points[a,0], points[b,0] ],
                       [ points[a,1], points[b,1] ], "k.-")
    
    plt.draw()
    plt.pause(0.001)

print()
plt.show()
