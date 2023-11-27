#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

R = 8.314
def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

G1 = lambda T, x: (1-x)*(-7000+120*T) + x*(-6000+130*T) + R*T*( xlnx(1-x) + xlnx(x) ) + (1-x)*x*(35000-10*T -4000*(1-x-x))
G2 = lambda T, x: (1-x)*(-5000+120*T) + x*(-8000+130*T) + R*T*( xlnx(1-x) + xlnx(x) ) + (1-x)*x*(30000-10*T)
G3 = lambda T, x: (1-x)*(+4000+110*T) + x*(+5000+120*T) + R*T*( xlnx(1-x) + xlnx(x) ) + (1-x)*x*(20000-5*T)

phases = [G1, G2, G3]

N = 250
x = np.linspace(0, 1, N)

points = np.empty( (len(phases)*N, 2) )
phase_id = np.empty(points.shape[0], dtype=int)

fig, ax = plt.subplots(1, 2, figsize=(14,6))
ax[1].grid()
ax[1].set_ylim(500, 1500)
ax[1].set_xlim(0, 1)

for T in range(500, 1500, 10):
    print(".", end="", flush=True)

    ax[0].cla()
    ax[0].grid()
    ax[0].set_xlim(0, 1)
    for G in phases:
        ax[0].plot(x, G(T, x), "-")

    o = 0
    for i,G in enumerate(phases):
        points[o:o+N,0] = x
        points[o:o+N,1] = G(T, x)
        phase_id[o:o+N] = i
        o += N

    hull = ConvexHull(points)
    mask = hull.equations[:,1] < 0 
    for a, b in hull.simplices[mask]:
        if phase_id[a] != phase_id[b]:
            ax[1].plot([ points[a, 0], points[b, 0] ], [ T, T ], "k-")
            ax[1].plot([ points[a, 0] ], [ T ], f"C{phase_id[a]}.")
            ax[1].plot([ points[b, 0] ], [ T ], f"C{phase_id[b]}.")
            ax[0].plot([ points[a,0], points[b,0] ],
                       [ points[a,1], points[b,1] ], "k.-")
    
    plt.draw()
    plt.pause(0.001)

print()
plt.show()
