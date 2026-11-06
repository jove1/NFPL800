#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from interactive_figures import InteractiveFigureSlope, InteractiveFigurePos

np.seterr(all='print')

R = 8.314

def xlnx(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x, x*np.log(x), 0)

class Phase2SL:
    def __init__(self, aa, ab, ba, bb):
        self.aa = aa
        self.ab = ab
        self.ba = ba
        self.bb = bb
    
    def X(self, y):
        y1a, y1b, y2a, y2b = y
        return np.array([y1a + y2a, 
                         y1b + y2b])/2

    def G(self, y):
        y1a, y1b, y2a, y2b = y
        return (
                y1a * y2a * self.aa + 
                y1a * y2b * self.ab + 
                y1b * y2a * self.ba +
                y1b * y2b * self.bb +
                RT*( xlnx(y1a) + xlnx(y1b) + xlnx(y2a) + xlnx(y2b) )
        )

    def dGdy(self, y):
        y1a, y1b, y2a, y2b = y
        return np.array([
            y2a * self.aa + y2b * self.ab + RT*( 1 + np.log(y1a) ),
            y2a * self.ba + y2b * self.bb + RT*( 1 + np.log(y1b) ),
            y1a * self.aa + y1b * self.ba + RT*( 1 + np.log(y2a) ),
            y1a * self.ab + y1b * self.bb + RT*( 1 + np.log(y2b) ),
        ])
    
    def dXdy(self, y):
        return np.array([[1, 0, 1, 0], 
                         [0, 1, 0, 1]])/2
        

RT = R * 100
#RT = R * 200
#RT = R * 200
aa, bb = 10000, 5000
ab = ba = (aa+bb)/2 - 5000
#ab = ba = (aa+bb)/2 + 5000
phase = Phase2SL(aa, ab, ba, bb) 


x = np.linspace(0, 1, 100)
y1b, y2b = np.meshgrid(x, x)
y = [1-y1b, y1b, 1-y2b, y2b]
xx = phase.X(y)
gg = phase.G(y)
print(xx.shape)
print(gg.shape)


from scipy import optimize

import time
started = time.time()
g = []
for x1 in x:
#PROBLEM
    # Minimize Gibbs energy function at each required composition, 
    # don't forget the appropriate constraints
    sol = optimize.minimize(lambda y: (y[0]-0.4)**2 +(y[2]-0.2)**2,
        #x0=[0.5, 0.5, 0.5, 0.5],
        x0=[0.4, 0.6, 0.6, 0.4],
        bounds=[(0, 1)]*4,
        constraints = [
           dict(type="eq", fun=lambda y: y[1] - 0.7),
           dict(type="eq", fun=lambda y: y[3] - 0.3) 
        ]
    )
#END
    print(sol.x, sol.message)
    g.append( sol.fun )
print("Minimisation done in", time.time()-started)

plt.figure()
plt.title("2SL minimize")
plt.plot(xx[1].ravel(), gg.ravel(), ".", ms=0.5)
plt.plot(x, g, ".-")
plt.xlim(0,1)


from scipy.spatial import ConvexHull
points = np.column_stack([ xx[1].ravel(), gg.ravel() ])

started = time.time()
hull = ConvexHull(points)
print("Convex hull in", time.time()-started)

mask = hull.equations[:,1] < 0
#PROBLEM
# Select correct lines for lower part of convex hull
lines = [points[::131]]
#END

plt.figure()
plt.title("2SL convex hull")
plt.plot(xx[1].ravel(), gg.ravel(), ".", ms=0.5)
from matplotlib.collections import LineCollection
plt.gca().add_collection( LineCollection(lines, color="C1") )
plt.xlim(0,1)



def calc_slope(slope):
    from scipy import optimize
#PROBLEM
    # Find point with given slope of phase.G as a minimum of a function,
    # don't forget to include correct constraints
    sol = optimize.minimize(lambda y: (y[0]-0.4)**2 + (y[2]-0.4)**2,
        x0=[0.4, 0.6, 0.4, 0.6],
        bounds=[(0, 1)]*4,
        constraints = [ 
           dict(type="eq", fun=lambda y: y[1] - 0.6),
           dict(type="eq", fun=lambda y: y[3] - 0.6) 
       ]
    )
#END
    print(sol.x, sol.message)
    y = sol.x

    g = phase.G(y)
    x = phase.X(y)
#PROBLEM
    # Calculate chemical potential at correct offset for phase.G
    mu = np.array([0, slope])
#END
    return [0, x[1], 1], [mu[0], g, mu[1]]

with InteractiveFigureSlope("2SL", calc_slope) as ax:
    ax.plot(xx[1].ravel(), gg.ravel(), ".", ms=0.5)
    ax.add_collection( LineCollection(lines, color="C2") )

def calc_pos(x1):
#PROBLEM
    g = phase.G([1, 0, 1-x1, x1])
    return [0, x1, 1], [ g, g, g ]
#END

with InteractiveFigurePos("2SL", calc_pos) as ax:
    ax.plot(xx[1].ravel(), gg.ravel(), ".", ms=0.5)
    ax.add_collection( LineCollection(lines, color="C2") )

plt.show()
