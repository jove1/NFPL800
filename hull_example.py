#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from matplotlib.collections import LineCollection

plt.rc("axes", grid=True)

# Generate random points inside a circle of radius 10.
points = np.random.uniform(-10, 10, size=(200, 2))
points = points[np.linalg.norm(points, axis=1) < 10]

# Compute the convex hull of the point cloud.
hull = ConvexHull(points)

# Plot the hull edges together with their outward normal vectors.
plt.figure()
plt.gca().set_aspect("equal")
plt.plot(points[:, 0], points[:, 1], "C0o")

# Each simplex is an edge (given by two point indices). Its equation is
# nx*x + ny*y + offset = 0; (nx, ny) is the outward normal.
for (a, b), (nx, ny, offset) in zip(hull.simplices, hull.equations):
    ax, ay = points[a]
    bx, by = points[b]
    cx, cy = (points[a] + points[b]) / 2

    plt.plot([ax, bx], [ay, by], "C1.-")
    plt.plot([cx, cx + nx], [cy, cy + ny], "C2-")


# The same hull boundary can be drawn more compactly with a LineCollection.
plt.figure()
plt.gca().set_aspect("equal")
plt.plot(points[:, 0], points[:, 1], "C0o")

# Gather all edge endpoints; lines has shape (number_of_edges, 2, 2).
lines = points[hull.simplices]
lc = LineCollection(lines, color="C1")
plt.gca().add_collection(lc)

plt.show()
