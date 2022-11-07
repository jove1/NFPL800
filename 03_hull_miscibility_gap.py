#!/usr/bin/python3
from pylab import *
from scipy.spatial import ConvexHull

from dataCuAg import G_AgCu_fcc, G_AgCu_liq, T0

N = 150
x, dx = linspace(0, 1, N, retstep=True)

T = T0 + 600
#T = T0 + 900

G_fcc =  G_AgCu_fcc(T, x)
plot(x, G_fcc, label="fcc")

G_liq =  G_AgCu_liq(T, x)
plot(x, G_liq, label="liquid")


points = concatenate([ 
    stack([x, G_fcc], axis=1),
    stack([x, G_liq], axis=1),
    [ (0, 1e5), 
      (1, 1e5) ]
])

tags = np.array([1]*N + [2]*N + [0]*2)

hull = ConvexHull(points)

for a,b in hull.simplices:
    if not tags[a] or not tags[b]:
        continue

    elif tags[a] != tags[b]:
        plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "C2.-")

    elif abs(points[a,0] - points[b,0]) > 1.1*dx:
        plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "C3.-")
    else:
        plot([ points[a,0], points[b,0] ],
             [ points[a,1], points[b,1] ], "C4.-")
 
legend()
xlim(0,1)
grid()
show()

