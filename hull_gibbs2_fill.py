#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

R = 8.314
G1 = lambda T, xa, xb: xa*(-7000+120*T) + xb*(-6000+130*T) + R*T*( xlnx(xa) + xlnx(xb) ) + xa*xb*(35000-10*T-4000*(xa-xb))
G2 = lambda T, xa, xb: xa*(-5000+120*T) + xb*(-8000+130*T) + R*T*( xlnx(xa) + xlnx(xb) ) + xa*xb*(30000-10*T)
G3 = lambda T, xa, xb: xa*(+4000+110*T) + xb*(+5000+120*T) + R*T*( xlnx(xa) + xlnx(xb) ) + xa*xb*(20000-5*T)

phases = [G1, G2, G3]

T = 950
N = 250
x = np.linspace(0, 1, N)

plt.figure()
for G in phases:
    plt.plot(x, G(T, 1-x, x), "-")
plt.grid()

points = np.empty( (len(phases)*N, 2) )
phase_id = np.empty(points.shape[0], dtype=int)
o = 0
for G in phases:
    # FILL IN points
    points[:] = np.random.normal((0.5, 100000), (0.5,10000), size=points.shape)
    o += N

hull = ConvexHull(points)
mask = hull.equations[:,1] < 0 
for a, b in hull.simplices[mask]:
    plt.plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "k.-")

plt.show()
