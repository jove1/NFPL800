#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

plt.figure()

points = np.random.normal(size=(300, 2))
plt.plot(points[:,0], points[:,1], "o")

hull = ConvexHull(points)

for (a, b), (c,d,e) in zip(hull.simplices, hull.equations):
    l, = plt.plot([ points[a,0], points[b,0] ],
                  [ points[a,1], points[b,1] ])

    mid = (points[a]+points[b])/2
    
    plt.plot([ mid[0], mid[0] + c ],
             [ mid[1], mid[1] + d ], "--", c=l.get_color() )

plt.gca().set_aspect("equal")


plt.show()

