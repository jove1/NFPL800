#!/usr/bin/python3
from pylab import *
from scipy.spatial import ConvexHull

from dataCuAg import G_AgCu_fcc, G_AgCu_liq, T0

N = 250
Ntemp = 100
x, dx = linspace(0, 1, N, retstep=True)

points = np.empty((2*N+2,2))
points[:N,0] = x
points[N:2*N,0] = x
points[2*N:] = [ (0, 1e5), (1, 1e5) ]

tags = np.array([1]*N + [2]*N + [0]*2)

ax1 = figure().gca() if Ntemp <= 10 else None
ax2 = figure().gca()

for T in linspace(20, 1200, Ntemp):
    print(".", end="", flush=True)
    points[:N,   1] =  G_AgCu_fcc(T0+T, 1-x)
    points[N:2*N,1] =  G_AgCu_liq(T0+T, 1-x)
   
    if ax1:
       ax1.plot( points[:N,0], points[:N,1], "C0-") 
       ax1.plot( points[N:2*N,0], points[N:2*N,1], "C1-")

    hull = ConvexHull(points)
    t = tags[hull.simplices]
    
    # phase equilibria
    mask = (t[:,0] != t[:,1]) & (t[:,0] != 0) & (t[:,1] != 0)
    for a,b in hull.simplices[mask]:
        ax2.plot([ points[a,0], points[b,0] ], [T, T], "C2.-")
        if ax1:
            ax1.plot([ points[a,0], points[b,0] ],
                     [ points[a,1], points[b,1] ], "C2.-")
    
    # miscibility gap in one phase
    mask = (t[:,0] == t[:,1]) & (t[:,0] != 0) & (t[:,1] != 0) & \
           (abs(points[hull.simplices[:,0],0] - \
                points[hull.simplices[:,1],0]) > 1.1*dx)
    for a,b in hull.simplices[mask]:
        ax2.plot([ points[a,0], points[b,0] ], [T, T], "C3.-")
        if ax1:
            ax1.plot([ points[a,0], points[b,0] ],
                     [ points[a,1], points[b,1] ], "C3.-")
print()
ax2.set_xlim(0,1)
ax2.grid(True)
show()

