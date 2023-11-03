#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

R = 8.314
N = 150
T = 350

x = np.linspace(0, 1, N)
xa, xb = 1-x, x
G = R*T*( xlnx(xa) + xlnx(xb) ) + xa*xb*(8000 - 2000*(xa-xb))

plt.figure()
plt.plot(x, G, ".")
plt.grid()

points = np.empty((N,2))
points[:,0] = np.random.normal(size=N) # FILL IN
points[:,1] = np.random.normal(0,100,size=N) # FILL IN

hull = ConvexHull(points)
for a, b in hull.simplices:
    plt.plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "-")


plt.show()

