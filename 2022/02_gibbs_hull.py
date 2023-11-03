#!/usr/bin/env python3
from pylab import *
from scipy.spatial import ConvexHull


def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return where(x, x*log(x), 0)

T0 = 273.15
R = 8.314
def Gmix(x, T):
    return T*R*( xlnx(x) + xlnx(1-x) )

N = 150
x = linspace(0, 1, N)


Ga =  Gmix(x, T0+20)  + 3000*x*(1-x)*(1-x-x)
Gb =  Gmix(x, T0+20)  - 4000*x*(1-x)*(1-x-x)
plot(x, Ga)
plot(x, Gb)

points = concatenate([ 
    stack([x, Ga], axis=1),
    stack([x, Gb], axis=1),
    [ (0, 1e5), 
      (1, 1e5) ]
])

tags = np.array([1]*N + [2]*N + [0]*2)

hull = ConvexHull(points)

for a,b in hull.simplices:
    if not tags[a] or not tags[b]:
        continue
    if tags[a] != tags[b]:
        plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "C2o-")
    elif tags[a] == tags[b] == 1:
        plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "C3.-")
    elif tags[a] == tags[b] == 2:
        plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "C4.-")

show()

