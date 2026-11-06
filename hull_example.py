#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

points = np.random.uniform(-10, 10, size=(200,2))
points = points[ np.linalg.norm(points, axis=1)<10 ]

hull = ConvexHull(points)

plt.figure()
plt.gca().set_aspect("equal")
plt.plot(points[:,0], points[:,1], "C0o")
for (a, b), (nx, ny, o) in zip(hull.simplices, hull.equations):
    ax, ay = points[a]
    bx, by = points[b]
    cx, cy = (points[a] + points[b])/2

    plt.plot([ax, bx], [ay, by], "C1.-")
    plt.plot([cx, cx+nx], [cy, cy+ny], "C2-")


plt.figure()
plt.gca().set_aspect("equal")
plt.plot(points[:,0], points[:,1], "C0o")
lines = points[hull.simplices]

from matplotlib.collections import LineCollection
lc = LineCollection(lines, color="C1")
plt.gca().add_collection(lc)


plt.show()

