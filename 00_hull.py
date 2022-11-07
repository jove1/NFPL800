#!/usr/bin/python3
from pylab import *
from scipy.spatial import ConvexHull

points = random((10,2))
hull = ConvexHull(points)

plot(points[:,0], points[:,1], "o")
for a,b in hull.simplices:
    plot([ points[a,0], points[b,0] ],
         [ points[a,1], points[b,1] ], "-")

show()

