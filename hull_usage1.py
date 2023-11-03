#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

plt.figure()

points = np.random.normal(size=(300, 2))
plt.plot(points[:,0], points[:,1], "o")

hull = ConvexHull(points)

# explicit:
for a, b in hull.simplices:
    plt.plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ])

# same as:
plt.plot(*points[hull.simplices].T)


plt.show()

